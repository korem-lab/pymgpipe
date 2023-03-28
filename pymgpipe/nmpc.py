import os
import tqdm
import logging
import pandas as pd
import optlang
import time
from collections import namedtuple
from pathlib import Path
from .fva import regularFVA, veryFastFVA
from .utils import load_dataframe, load_model, set_objective, Constants
from .io import suppress_stdout

def compute_nmpcs(
    samples,
    out_dir="./",
    out_file="nmpcs.csv",
    objective_out_file="community_objectives.csv",
    fluxes_out_file="all_fluxes.csv",
    reactions=None,
    regex=None,
    ex_only=True,
    solver="gurobi",
    threads=int(os.cpu_count() / 2),
    parallel=True,
    diet_fecal_compartments=True,
    force=False,
    threshold=1e-5,
    write_to_file=True,
    fva_type="regular",
    obj_optimality=100,
    scaling=0,
):
    """Compute NMPCs as well as associated reaction metrics on specified list (or directory) of samples

    This function is computes NMPCs at the level of the optlang LP problem to leverage speed of low-level representation.
    Available FVA types are `regular` and `fast`. `Fast` FVA is significantly faster, but requires both a CPLEX license and full install of the VFFVA package.
    Additionally, all problems need to be in .mps format if using `fast` FVA type. For more information, see VFFVA package- https://github.com/marouenbg/VFFVA.

    NMPCs are calculated as the max fluxes secreted by the fecal compartment (positive values only) summed with the min flux uptaken by the diet compartment (negative values only).
    In cases where models are not build with diet/fecal compartments, NMPCs are simply calculated as the sum of the max and min fluxes.

    Args:
        samples (list | str): List of samples or directory containing samples
        out_dir (str): Directory to output results
        out_file (str): Name of file containing final NMPCs
        objective_out_file (str): Name of file containing community objectives
        fluxes_out_file (str): Name of file containing fluxes
        reactions (list): List of reactions to run NMPCs on
        regex (str): Regex match for list of reactions to run NMPCs on
        diet_fecal_compartments (bool): Whether or not models are built with diet/fecal compartmentalization 
        ex_only (bool): Compute NMPCs on exchange reactions only
        fva_type (str): FVA type used to compute NMPCs, allowed values are `fast` and `regular`
        solver (str): LP solver used to compute NMPCs, allowed values are `gurobi` and `cplex`
        obj_optimality (float): Percent of optimal objective value constrained during NMPC computation
        threshold (float): Fluxes below threshold will be set to 0
        write_to_file (bool): Write results to file

    Notes:
        If computation is cut short prematurely, this function will pick up where it left off based on which samples are already present in `out_file`. 

    """
    assert fva_type == "regular" or fva_type == "fast", (
        "FVA type must be either `regular` or `fast`! Received %s" % fva_type
    )

    start = time.time()
    out_dir = out_dir + "/" if out_dir[-1] != "/" else out_dir
    Path(out_dir).mkdir(exist_ok=True)

    out_file = out_dir + out_file
    objective_out_file = out_dir + objective_out_file
    fluxes_out_file = out_dir + fluxes_out_file

    nmpcs = (
        pd.DataFrame()
        if force or not write_to_file
        else load_dataframe(out_file, return_empty=True)
    )
    all_fluxes = (
        pd.DataFrame()
        if force or not write_to_file
        else load_dataframe(fluxes_out_file, return_empty=True)
    )
    obj_values = (
        pd.DataFrame()
        if force or not write_to_file
        else load_dataframe(objective_out_file, return_empty=True)
    )
    obj_values["communityBiomass"] = (
        None if obj_values.empty else obj_values[obj_values.columns[0]]
    )

    try:
        models = [load_model(samples)]

        # Skip models that already exist
        if models[0].name in list(nmpcs.columns) and not force:
            print("NMPCs for %s already exist in file!" % models[0].name)
            return
    except Exception:
        models = (
            samples
            if isinstance(samples, list)
            else [
                os.path.dirname(samples) + "/" + m
                for m in os.listdir(os.path.dirname(samples))
            ]
        )

        # Skip models that already exist
        models = [
            f
            for f in models
            if not isinstance(f, str)
            or os.path.basename(f).split(".")[0] not in list(nmpcs.columns)
        ]
    threads = os.cpu_count() - 1 if threads == -1 else threads

    print("Computing NMPCs on %s models..." % len(models))
    print("\n----------------Parameters----------------")
    print("Parallel- %s" % str(parallel).upper())
    print("Threads- %s" % str(threads).upper())
    print("Solver- %s" % solver.upper())
    print("Diet/fecal compartments- %s" % str(diet_fecal_compartments).upper())
    print("Exchanges only- %s" % str(ex_only).upper())
    print("Out file- %s" % str(out_file).upper())
    print("Force ovewrite- %s" % str(force).upper())
    print("------------------------------------------")

    for s in tqdm.tqdm(models, total=len(models)):
        if fva_type == "fast":
            solver = 'cplex'
            assert isinstance(
                s, str
            ), "For fast fva, `samples` param must be directory or list of model paths."
            model_path = s
        with suppress_stdout():
            m = load_model(path=s, solver=solver)
            if not isinstance(m, optlang.interface.Model):
                raise Exception("Expected optlang.Model, received %s" % type(m))
            if not force and m.name in list(nmpcs.columns):
                continue

        # Solve for objective first
        if "communityBiomass" not in m.variables:
            raise Exception("Could not find communityBiomass variable in model!")
        m.variables["communityBiomass"].set_bounds(0.4, 1)
        set_objective(m, m.variables["communityBiomass"], direction="max")

        # Now perform FVA under constrained objective value
        try:
            if fva_type == "regular":
                with suppress_stdout():
                    m.optimize()
                if m.status == "infeasible":
                    logging.warning("%s model is infeasible!" % m.name)
                    continue

                obj_val = round(m.objective.value, 5)
                obj_values.loc[m.name] = obj_val
                if "ObjectiveConstraint" in m.constraints:
                    m.remove(m.constraints["ObjectiveConstraint"])
                    m.update()
                obj_const = m.interface.Constraint(
                    expression=m.objective.expression,
                    lb=obj_val * (obj_optimality / 100),
                    ub=obj_val,
                    name="ObjectiveConstraint",
                )
                m.add(obj_const)
                m.update()
                res = regularFVA(
                    m,
                    reactions=reactions,
                    regex=regex,
                    ex_only=ex_only,
                    solver=solver,
                    threads=threads if parallel else 1,
                    parallel=parallel,
                    write_to_file=False,
                    threshold=threshold,
                )
            elif fva_type == "fast":
                res = veryFastFVA(
                    model=m,
                    path=model_path,
                    reactions=reactions,
                    regex=regex,
                    nCores=threads if parallel else 1,
                    nThreads=1,
                    optPerc=obj_optimality,
                    threshold=threshold,
                    scaling=scaling
                )
        except Exception as e:
            logging.warning(f"Cannot solve {m.name} model!\n{e}")
            continue
        if res is None:
            return
        res["sample_id"] = m.name
        all_fluxes = pd.concat([all_fluxes, res], axis=0)
        if diet_fecal_compartments:
            metabs = [
                m
                for m in res.index.str.split("[").str[0].drop_duplicates()
                if not m.startswith("Diet")
            ]
            df = {}
            for metab in metabs:
                fe = res.loc[metab + "[fe]"]["max"]
                d = res.loc["Diet_" + metab + "[d]"]["min"]

                df[metab.split("EX_")[1]] = d + fe
            nmpc = pd.DataFrame({m.name: df})
        else:
            nmpc = res["min"] + res["max"]
            nmpc.name = m.name

        nmpcs = pd.concat([nmpcs, nmpc], axis=1).fillna(0)
        if write_to_file:
            nmpcs.to_csv(out_file)
            obj_values.to_csv(objective_out_file)
            all_fluxes.to_csv(fluxes_out_file)

    res = namedtuple("res", "nmpc objectives fluxes")

    print("-------------------------------------------------------")
    print("Finished computing NMPCs!")
    print("Process took %s minutes to run..." % round((time.time() - start) / 60, 3))

    return res(nmpcs, obj_values, all_fluxes)