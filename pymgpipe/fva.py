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
    suppress_stdout
)

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

def regularFVA(
    model=None,
    reactions=None,
    regex=None,
    ex_only=True,
    solver="gurobi",
    threads=int(os.cpu_count() / 2),
    write_to_file=False,
    out_dir="fva/",
    parallel=True,
    threshold=1e-5,
):
    gc.enable()

    model = load_model(path=model, solver=solver)

    if reactions is None and regex is None and ex_only is True:
        regex = Constants.EX_REGEX

    reactions_to_run = [r.name for r in get_reactions(model, reactions, regex)]

    result_df = []
    if write_to_file:
        Path(out_dir).mkdir(exist_ok=True)
        out_file = out_dir + "%s.csv" % model.name

        if os.path.exists(out_file):
            result_df = pd.read_csv(out_file)
            result_df = result_df[~result_df["id"].isnull()][
                ~result_df.id.duplicated(keep="first")
            ]
            metabs_to_skip = list(result_df.id)
            print(
                "\nFound existing file, skipping %s metabolites..."
                % len(metabs_to_skip)
            )

            reactions_to_run = [r for r in reactions_to_run if r not in metabs_to_skip]
            result_df = result_df.to_dict("records")

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
    print(
        "Starting parallel FVA with %s chunks on %s threads"
        % (threads, len(split_reactions))
    )

    _func = partial(_optlang_worker, threshold)
    if parallel:
        p = Pool(processes=threads, initializer=partial(_pool_init, model))
        res = p.imap(_func, split_reactions)
    else:
        _pool_init(model)
        res = map(_func, split_reactions)

    for result in res:
        result_df = result_df + result

        out_df = pd.DataFrame.from_records(result_df, index="id")
        out_df.sort_index(inplace=True)
        if write_to_file:
            out_df.to_csv(out_file)

    try:
        p.close()
        p.join()
    except Exception:
        pass

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
        max_sol = solve_model(model=global_model, reactions=[m], flux_threshold=threshold).to_dict()[
            global_model.name
        ][m]

        global_model.objective = Objective(net, direction="min")
        min_sol = solve_model(model=global_model, reactions=[m], flux_threshold=threshold).to_dict()[
            global_model.name
        ][m]

        result.append({"id": m, "min": min_sol, "max": max_sol})
    return result

def veryFastFVA(
        nCores, 
        nThreads,
        path, 
        model=None, 
        ex_only=True,
        scaling=0, 
        memAff='none', 
        schedule='dynamic', 
        nChunk=50, 
        optPerc=90, 
        reactions=None,
        regex=None,
        threshold=1e-5):
    '''
    VFFVA performs Very Fast Flux Variability Analysis (VFFVA). VFFVA is a parallel implementation of FVA that
    allows dynamically assigning reactions to each worker depending on their computational load
    Guebila, Marouen Ben. "Dynamic load balancing enables large-scale flux variability analysis." bioRxiv (2018): 440701.

    USAGE:
    minFlux,maxFlux=VFFVA(nCores, nThreads, model, scaling, memAff, schedule, nChunk, optPerc, ex)

    :param nCores:   Number of non-shared memory cores/machines.
    :param nThreads: Number of shared memory threads in each core/machine.
    :param model:    .mps format: path to metabolic model in .mps format.
    :param scaling:  CPLEX parameter. It corresponds to SCAIND parameter (Default = 0).
                     -1: no scaling; 0: equilibration scaling; 1: more aggressive scaling.
                     more information here: https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/ScaInd.html.
    :param memAff:   'none', 'core', or 'socket'. (Default = 'none'). This an OpenMPI parameter, more
                     information here: https://www.open-mpi.org/faq/?category=tuning#using-paffinity-v1.4.
    :param schedule: 'dynamic', 'static', or 'guided'. (Default = 'dynamic'). This is an OpenMP parameter, more
                     information here: https://software.intel.com/en-us/articles/openmp-loop-scheduling
    :param nChunk:   Number of reactions in each chunk (Default = 50). This is an OpenMP parameter, more
                     information here: https://software.intel.com/en-us/articles/openmp-loop-scheduling
    :param optPerc:  Percentage of the optimal objective used in FVA. Integer between 0 and 100. For example, when set to 90
                     FVA will be computed on 90% of the optimal objective.
    :param ex:      0-based indices of reactions to optimize as a numpy array.  (Default, all reactions)
    :return:
           minFlux:          (n,1) vector of minimal flux values for each reaction.
           maxFlux:          (n,1) vector of maximal flux values for each reaction.
    '''
    status = os.system('mpirun --version')
    if status != 0:
        raise ValueError(['MPI and/or CPLEX nont installed, please follow the install guide'
               'or use the quick install script'])

    # Set schedule and chunk size parameters
    os.environ["OMP_SCHEDUELE"] = schedule+str(nChunk)

    model = load_model(path if model is None else model)
    with suppress_stdout():
        model.optimize()
    if model.status == "infeasible":
        raise Exception("%s model is infeasible!" % model.name)
        return
    
    var_dict = {i:v.name for i,v in enumerate(model.variables)}
    var_dict_inv = {v:k for k,v in var_dict.items()}

    ex_indices = [var_dict_inv[r.name] for r in get_reactions(model, reactions = reactions, regex = Constants.EX_REGEX if ex_only and regex is None else regex)]
    print('Running on %s reactions!'%len(ex_indices))
    
    # Set reactions to optimize
    rxns_file = 'rxns_%s.txt'%model.name
    with open(rxns_file, "w") as f:
        for num in ex_indices:
            f.write(str(num) + "\n")   
  
    assert vffva_config.path is not None, 'Please set value of `vffva_config.path` to location of VFFVA executable'
    status = os.system('mpirun -np ' + str(nCores) + ' --bind-to ' + str(memAff) + ' -x OMP_NUM_THREADS=' + str(nThreads) +
        f' {vffva_config.path}/lib/veryfastFVA ' + path + ' ' + str(optPerc) + ' ' + str(scaling) + ' ' + rxns_file)

    # Fetch results
    resultFile = path[:-4] + 'output.csv'
    if not os.path.exists(resultFile):
        raise Exception('Ran into issue when running VFFVA, could not find results file...')
    results = pd.read_csv(resultFile,header=0)

    os.system('rm '+resultFile)
    os.system('rm '+rxns_file)
    results.rename({i:var_dict[x] for i,x in enumerate(ex_indices) if i in results.index},axis=0,inplace=True)

    results.index.rename('id',inplace=True)
    if threshold is not None:
        results[abs(results)<threshold]=0
    
    results.sort_index(inplace=True)
    results.columns=['min','max']
    return results


def _pool_init(sample_model):
    sys.stdout = open(os.devnull, "w")

    global global_model
    global_model = sample_model

    global global_problem
    global_problem = sample_model.problem
