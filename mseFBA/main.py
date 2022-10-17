import os
import sys
import pandas as pd
from multiprocessing import Pool
from functools import partial
import tqdm
import gc
import numpy as np

from pymgpipe import load_model, solve_model
from pymgpipe.optlang_util import _get_reverse_id, get_reactions

from .metabolomics import *
from .utils import *

from optlang.interface import Objective
import time
import logging
from random import shuffle

class Constants:
    EX_REGEX = '^EX_.*_m$'

def run(
    metabolomics,
    samples='problems/',
    out_file=None,
    dataset_dir='./',
    transformation=transformations.none,
    ex_only=True,
    zero_unmapped_metabolites=False,
    threads=int(os.cpu_count()/2),
    solver='gurobi',
    verbosity=0,
    presolve=True,
    threshold=1e-5,
    parallel=True,
    metabolites=None,
):
    start_time = time.time()
    gc.disable()
 
    try:
        model_files = samples if isinstance(samples,list) else [dataset_dir+samples+m for m in os.listdir(dataset_dir+samples)]
    except:
        raise Exception('Please pass in a valid model directory or an explicit list of sample paths using the \'problems\' parameter')

    metabolomics_df = transform_metabolomics(metabolomics,transformation)
    unmatched_metabolomics = [f for f in model_files if f.split('/')[-1].split('.')[0] not in list(metabolomics_df.columns)]
    if len(unmatched_metabolomics) > 0:
        print('%s samples dont have associated columns in metabolomics file-\n'%len(unmatched_metabolomics))
        print([f.split('/')[-1].split('.')[0] for f in unmatched_metabolomics])
        model_files = [f for f in model_files if f not in unmatched_metabolomics]
        
    final = load_dataframe(out_file,return_empty=True)
    if len(final.columns)>0:
        print('Skipping %s samples that are already finished!'%len(final.columns))
        model_files = [f for f in model_files if f.split('/')[-1].split('.')[0] not in list(final.columns)]
        shuffle(model_files)
        
    if len(model_files) == 0 and final.empty:
        logging.warning('Something went wrong! No existing solution or samples left to run.')
        return (None, None)
    elif len(model_files) == 0:
        print('Finished mseFBA, no samples left to run!')
        print('The results are in...\n')
        res = evaluate_results(final,metabolomics_df)
        return (final,res)

    threads = os.cpu_count()-1 if threads == -1 else threads
    threads = min(threads,len(model_files))

    print('\n----------------Parameters----------------')
    print('Parallel- %s'%str(parallel).upper())
    print('Solver- %s'%solver.upper())
    print('Presolve- %s'%str(presolve).upper())
    print('Verbosity- %s'%verbosity)
    print('Zero unmapped metabolites- %s'%str(zero_unmapped_metabolites).upper())
    print('------------------------------------------')
        
    _func = partial(
        _mseFBA_worker,
        ex_only,
        zero_unmapped_metabolites,
        solver,
        verbosity,
        presolve,
        threshold,
        metabolites
    )
    
    if parallel:
        print('Running mseFBA on %s samples in parallel using %s threads...\n'%(len(model_files),threads))
        p = Pool(processes=threads,initializer=partial(_pool_init,metabolomics_df))
        res = p.imap(_func, model_files)
    else:
        print('Running mseFBA on %s samples in series...\n'%len(model_files))
        _pool_init(metabolomics_df)
        res = map(_func, model_files)
    
    infeasible = []
    feasible = []
    objective_vals = {}
    for r in tqdm.tqdm(res,total=len(model_files)):
        sample_file, solution, obj = r
        if solution is not None:
            final = pd.concat([final,solution],axis=1)
            final.sort_index(inplace=True)
            if out_file is not None:
                final.to_csv(out_file)
            feasible.append(sample_file)
            objective_vals[solution.columns[0]]=obj
        else:
            infeasible.append(sample_file)
    try:
        p.close()
        p.join()
    except:
        pass

    sys.stdout = sys.__stdout__
    print('\n--Finished mseFBA in %s seconds! Solved %s samples and saved to %s--'%((int(time.time()-start_time)),len(feasible),out_file))
    if len(infeasible) > 0:
        print('Unable to solve %s models-\n'%len(infeasible))
        print(infeasible)

    print('The results are in...\n')
    res = evaluate_results(final,metabolomics_df)
    objectives = pd.DataFrame({'objective':objective_vals})
    return (final,res,objectives)

def _mseFBA_worker(ex_only, zero_unmapped_metabolites, solver, verbosity, presolve, threshold, metabolites, model_file):
    global solution_path, metabolomics_global
    model = load_model(model_file,solver)
    if model.name not in metabolomics_global.columns:
        raise Exception('Could not find %s in metabolomics- %s'%(model_file,list(metabolomics_global.columns)))

    all_metabolites = metabolomics_global[model.name].dropna().to_dict()
    metabs_to_map = {k:v for k,v in all_metabolites.items() if k in metabolites} if metabolites is not None else all_metabolites
        
    if zero_unmapped_metabolites:
        ex_reactions = get_reactions(model, regex = Constants.EX_REGEX)
        unmapped_metabs = [k.name for k in ex_reactions if k.name not in all_metabolites]
        metabs_to_map.update({k:0 for k in unmapped_metabs})
        
    add_correlation_objective(model,metabs_to_map)

    solution = None
    obj_val = None
    try:
        solution = solve_model(
            model=model,
            regex = Constants.EX_REGEX if ex_only else None,
            verbosity=verbosity,
            presolve=presolve,
            method='barrier',
            flux_threshold=threshold
        )
        obj_val = get_objective_value(model)
    except:
        pass

    del model
    # gc.collect()
    
    return (model_file, solution, obj_val)
 
def add_correlation_objective(model, flux_map):
    obj_expression = get_mse_expression(model, flux_map)

    if obj_expression is None:
        logging.warning('No metabolites in correlation objective, returning model as-is.')
        return
    try:
        model.objective = Objective(obj_expression,direction="min")
        model.update()
    except Exception as e:
        raise Exception('Failed to add mseFBA objective to model- %s'%e)


def get_mse_expression(model, flux_map, multi_sample=False):
    obj_expression = None

    flux_map = {k:v for k,v in flux_map.items() if k in model.variables}
    for f_id, flux in flux_map.items():
        forward_var = model.variables[f_id]
        reverse_var = model.variables[_get_reverse_id(f_id, multi_sample)]
        net = forward_var-reverse_var

        squared_diff=(net-flux)**2
        obj_expression = squared_diff if obj_expression is None else obj_expression + squared_diff
    return obj_expression

def get_variance_expression(model, ids, multi_sample=False):
    vrs = [model.variables[f_id]-model.variables[_get_reverse_id(f_id, multi_sample)] for f_id in ids if f_id in model.variables]
    if len(vrs) <= 1:
        return None
    mean = None
    for v in vrs:
        mean = v if mean is None else mean + v    
    mean = mean/len(vrs)

    obj_expr = np.sum([(v-mean)*(v-mean) for v in vrs])
    return obj_expr

def _pool_init(m_df):
    sys.stdout = open(os.devnull, 'w')  

    global metabolomics_global
    metabolomics_global = m_df
