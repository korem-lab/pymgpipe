import os
import gc
import sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
from functools import partial
from optlang.interface import Objective
from pathlib import Path
from .utils import (
    load_model,
    Constants,
    get_reactions,
    get_reverse_id,
    solve_model,
    load_dataframe,
)
from .io import suppress_stdout
from .logger import logger
from enum import Enum

class FVA_TYPE(Enum):
    REGULAR = 1
    FAST = 2

class VFFVA(object):
    def __init__(self):
        self._path = None

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        self._path = value

vffva_config = VFFVA()

def fva(
    model=None,
    fva_type=FVA_TYPE.REGULAR,
    reactions=None,
    regex=None,
    ex_only=True,
    solver="gurobi",
    threads=int(os.cpu_count() / 2),
    write_to_file=False,
    out_file=None,
    parallel=True,
    threshold=1e-5,
    scaling=0,
    mem_aff="none",
    schedule="dynamic",
    objective_percent=None,
    force=False,
):
    gc.enable()

    if fva_type == FVA_TYPE.FAST:
        assert isinstance(
            model, str
        ), "For fast FVA, `model` needs to be passed in as path to .mps file"
        assert os.path.isdir(vffva_config.path), 'Could not find VFFVA folder at %s'%vffva_config.path
        path = model

    model = load_model(path=model, solver=solver)
    with suppress_stdout():
        model.optimize()
    if model.status == "infeasible":
        raise Exception("%s model is infeasible!" % model.name)

    if reactions is None and regex is None and ex_only is True:
        regex = Constants.EX_REGEX

    reactions_to_run = [r.name for r in get_reactions(model, reactions, regex)]

    out_df = pd.DataFrame()
    if write_to_file:
        out_file = out_file if out_file is not None else "%s_fva.csv" % model.name
        out_df = load_dataframe(out_file, return_empty=True) if not force else pd.DataFrame()

        metabs_to_skip = list(out_df.index)
        if len(metabs_to_skip) > 0:
            print("\nFound existing file, skipping %s metabolites..." % len(metabs_to_skip))
            reactions_to_run = [r for r in reactions_to_run if r not in metabs_to_skip]

    if len(reactions_to_run) == 0:
        print("---Finished! No reactions left to run---")
        return

    if parallel is False:
        threads = 1
        parallel = False
    else:
        threads = (
            os.cpu_count() if threads == -1 or threads > os.cpu_count() else threads
        )
        threads = min(threads, len(reactions_to_run))

        parallel = False if threads <= 1 else parallel

    split_reactions = np.array_split(reactions_to_run, threads)
    logger.info(
        "Starting parallel FVA on %s reactions using %s chunks on %s threads...\n"
        % (len(reactions_to_run), threads, len(split_reactions))
    )

    if fva_type == FVA_TYPE.REGULAR:
        # VFFVA does this automatically
        if objective_percent is not None:
            print(f'Constraining objective value to {objective_percent}% of optimal value...')
            
            prev_bounds = {}
            primals = {v:v.primal for v in model.objective.expression.free_symbols}
            for v, value in primals.items():
                prev_bounds[v.name]=(v.lb,v.ub)
                v.lb = float(value) * (objective_percent/100)
                v.ub = float(value)
                print(v)
            
            model.update()

        _func = partial(_optlang_worker, threshold)
        if parallel:
            p = Pool(processes=threads, initializer=partial(_pool_init, model))
            res = p.imap(_func, split_reactions)
        else:
            _pool_init(model)
            res = map(_func, split_reactions)

        for result in res:
            result_df = pd.DataFrame.from_records(result, index="id")
            out_df = pd.concat([out_df, result_df], axis=0)
            out_df.sort_index(inplace=True)

            if write_to_file:
                out_df.to_csv(out_file)
        try:
            p.close()
            p.join()
        except Exception:
            pass

        # reset
        if objective_percent is not None:
            for v, (lower, upper) in prev_bounds.items():
                model.variables[v].set_bounds(lower, upper)

            model.update()

    elif fva_type == FVA_TYPE.FAST:
        status = os.system("mpirun --version")
        if status != 0:
            raise ValueError(
                [
                    "MPI and/or CPLEX nont installed, please follow the install guide"
                    "or use the quick install script"
                ]
            )

        # Set schedule and chunk size parameters
        os.environ["OMP_SCHEDUELE"] = schedule + str(len(split_reactions[0]))

        objective_percent = objective_percent if objective_percent is not None else -1

        var_dict = {i: v.name for i, v in enumerate(model.variables)}
        var_dict_inv = {v: k for k, v in var_dict.items()}

        ex_indices = [var_dict_inv[r] for r in reactions_to_run]

        # Set reactions to optimize
        rxns_file = "rxns_%s.txt" % model.name
        with open(rxns_file, "w") as f:
            for num in ex_indices:
                f.write(str(num) + "\n")

        assert vffva_config.path is not None, "Please set value of `vffva_config.path` to location of VFFVA executable"
        
        try:
            status = os.system( "mpirun -np " + str(threads) + " --bind-to " + str(mem_aff) + " -x OMP_NUM_THREADS=" + str(1) + f" {vffva_config.path}/lib/veryfastFVA " + path + " " + str(objective_percent) + " " + str(scaling) + " " + rxns_file )
        except:
            raise Exception(
                "Ran into issue when submitting VFFVA mpirun, please check installation- %s"
                % status
            )

        # Fetch results
        resultFile = path[:-4] + "output.csv"
        if not os.path.exists(resultFile):
            raise Exception(
                "Ran into issue when running VFFVA, could not find results file..."
            )
        result_df = pd.read_csv(resultFile, header=0)

        os.system("rm " + resultFile)
        os.system("rm " + rxns_file)
        result_df.rename(
            {i: var_dict[x] for i, x in enumerate(ex_indices) if i in result_df.index},
            axis=0,
            inplace=True,
        )
        result_df.index.rename("id", inplace=True)
        result_df.columns = ['min', 'max']
        
        out_df = pd.concat([out_df, result_df], axis=0)
    else:
        raise Exception('`fva_type` must be equal to either FVA_TYPE.REGULAR or FVA_TYPE.FAST! Received %s'%fva_type)

    if threshold is not None:
        out_df[abs(out_df) < threshold] = 0

    out_df.sort_index(inplace=True)
    out_df.columns = ['min', 'max']

    if write_to_file:
        out_df.to_csv(out_file)

    return out_df


# works for both cplex and gurobi
def _optlang_worker(threshold, metabolites):
    global global_model

    result = []
    for m in metabolites:
        net = global_model.variables[m]

        reverse_id = get_reverse_id(m)
        if reverse_id in global_model.variables:
            net -= global_model.variables[reverse_id]

        global_model.objective = Objective(net, direction="max")
        max_sol = solve_model(model=global_model, reactions=[m]).to_dict()[
            global_model.name
        ][m]

        global_model.objective = Objective(net, direction="min")
        min_sol = solve_model(model=global_model, reactions=[m]).to_dict()[
            global_model.name
        ][m]

        result.append({"id": m, "min": min_sol, "max": max_sol})
    return result

def _pool_init(sample_model):
    sys.stdout = open(os.devnull, "w")

    global global_model
    global_model = sample_model

    global global_problem
    global_problem = sample_model.problem
