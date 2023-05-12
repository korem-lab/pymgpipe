import os
import tqdm
import logging
import pandas as pd
import optlang
import time
from collections import namedtuple
from pathlib import Path
from .fva import FVA_TYPE, fva
from .utils import load_dataframe, load_model, set_objective, Constants
from .io import suppress_stdout
import cobra
import optlang

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
    fva_type=FVA_TYPE.REGULAR,
    objective_percent=100,
    scaling=0,
    mem_aff="none",
    schedule="dynamic",
    signed=False
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

    if isinstance(samples, str) and os.path.isdir(samples):
        models = [
            os.path.dirname(samples) + "/" + m
            for m in os.listdir(os.path.dirname(samples))
        ]
    elif isinstance(samples, list):
        models = samples 
    else:
        models = [samples]
  
    models = [
        f for f in models if force
        or (
            f.split("/")[-1].split(".")[0] not in list(nmpcs.columns)
            if isinstance(f, str) else f.name not in list(nmpcs.columns)
        )
    ]
    print("Computing NMPCs on %s models using %s..." % (len(models), str(fva_type)))

    for m in tqdm.tqdm(models, total=len(models)):
        m_name = m.split("/")[-1].split(".")[0] if isinstance(m, str) else m.name
        try:
            res = fva(
                m,
                solver=solver,
                fva_type=fva_type,
                reactions=reactions,
                regex=regex,
                ex_only=ex_only,
                threads=threads,
                parallel=parallel,
                write_to_file=False,
                threshold=threshold,
                objective_percent=objective_percent,
                scaling=scaling,
                mem_aff=mem_aff,
                schedule=schedule,
            )
        except Exception as e:
            logging.warning(f"Cannot solve {m_name} model!\n{e}")
            continue
        if res is None:
            return
        res["sample_id"] = m_name
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
            nmpc = pd.DataFrame({m_name: df})
        else:
            nmpc = res["min"] + res["max"]
            nmpc.name = m_name

        nmpcs = pd.concat([nmpcs, nmpc], axis=1).fillna(0)
        if not signed:
            nmpcs = abs(nmpcs)

        if write_to_file:
            nmpcs.to_csv(out_file)
            obj_values.to_csv(objective_out_file)
            all_fluxes.to_csv(fluxes_out_file)

    res = namedtuple("res", "nmpc objectives fluxes")

    print("-------------------------------------------------------")
    print("Finished computing NMPCs!")
    print("Process took %s minutes to run..." % round((time.time() - start) / 60, 3))

    return res(nmpcs, obj_values, all_fluxes)
