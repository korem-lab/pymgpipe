import os

import pandas as pd
from gurobipy import *

from multiprocessing import Pool
from functools import partial
import tqdm
import gc

from pymgpipe import InfeasibleModelException, load_model, solve_model
from pymgpipe.optlang_util import _get_exchange_reactions, _get_reverse_id

import numpy as np
import sys
from pathlib import Path

def run(
    metabolomics,
    files=[],
    file_dir='problems/',
    ex_only=True,
    zero_unmapped_metabolites=False,
    threads=int(os.cpu_count()/2),
    solver='gurobi',
    verbosity=0,
    presolve=True,
    threshold=1e-5,
    scale=True,
    map_labels=True,
    out_dir='mse/',
    parallelize=True,
    out_file=None
):
    gc.enable()

    Path(out_dir).mkdir(exist_ok=True)
        
    try:
        model_files = files if files else [file_dir+m for m in os.listdir(file_dir)]
    except:
        raise Exception('Please pass in a valid model directory using the \'dir\' parameter or an explicit set of files using the \'files\' parameter')

    finished = [f.split('.csv')[0] for f in os.listdir(out_dir)]
    if len(finished)>0:
        print('Skipping %s samples that are already finished!'%len(finished))
        model_files = [f for f in model_files if f.split('/')[-1].split('.mps')[0] not in finished]

    if len(model_files) == 0:
        print('Finished mseFBA, no samples left to run!')
        return

    threads = os.cpu_count() if threads == -1 else threads
    threads = min(threads,len(model_files))

    metabolomics_df = _get_metabolomics(metabolomics, scale, map_labels)

    print('Parallel- %s'%str(parallelize).upper())
    print('Solver- %s'%solver.upper())
    print('Presolve- %s'%str(presolve).upper())
    print('Verbosity- %s'%verbosity)
    print('Zero unmapped metabolites- %s'%str(zero_unmapped_metabolites).upper())
    print('Scale metabolomics- %s'%str(scale).upper())
    print('Map sample IDs- %s'%str(map_labels).upper())
    print('------------------------------------------')
        
    _func = partial(
        _mseFBA_worker,
        ex_only,
        zero_unmapped_metabolites,
        solver,
        verbosity,
        presolve,
        threshold,
        out_dir
    )
    if parallelize:
        print('Running mseFBA on %s samples in parallel using %s threads...\n'%(len(model_files),threads))
        p = Pool(processes=threads,initializer=partial(_pool_init,metabolomics_df))
        res = list(p.imap(_func, model_files))
        p.close()
        p.join()
    else:
        print('Running mseFBA on %s samples in series...\n'%len(model_files))
        _pool_init(metabolomics_df)
        res = list(map(_func,model_files))

    feasible_models = list(filter(lambda sample: sample[1] == True, res))
    infeasible_models = list(filter(lambda sample: sample[1] == False, res))

    print('Finished mseFBA! Solved %s samples'%len(feasible_models))
    if len(infeasible_models) > 0:
        print('Some models were infeasible and could not be solved-\n')
        print(infeasible_models)

    if out_file is not None:
        combine_mse_solutions(mse_dir=out_dir,out_file=out_file)

def _mseFBA_worker(ex_only, zero_unmapped_metabolites, solver, verbosity, presolve, threshold, out_dir, model_file):
    model = load_model(model_file,solver)

    metab_map = metabolomics_global[model.name].dropna().to_dict()
    if zero_unmapped_metabolites:
        ex_reactions = _get_exchange_reactions(model)
        metab_map.update({k:0 for k in ex_reactions if k not in metab_map})
        
    _add_correlation_objective(model,metab_map)

    solution = None
    try:
        solution = solve_model(
            model=model,
            ex_only=ex_only,
            verbosity=verbosity,
            presolve=presolve,
            method='barrier',
            flux_threshold=threshold
        )
    except InfeasibleModelException:
        raise Exception('UH OH')

    solved = solution is not None
    if solved:
        solution.sort_index(axis=0,inplace=True)
        solution.to_csv(out_dir+'%s.csv'%model.name)

    del model
    gc.collect()
    
    return (model_file, solved)
    
def _get_metabolomics(metabolomics,scale=True,map_labels=True):
    metabolomics_df = _load_dataframe(metabolomics)
    
    if scale:
        metabolomics_df = scale_metabolomics(metabolomics_df)

    if map_labels:
        conversion = _load_dataframe('sample_label_conversion.csv').conversion.to_dict()
        metabolomics_df.rename(conversion,axis='columns',inplace=True)

    return metabolomics_df

def _pool_init(m_df):
    sys.stdout = open(os.devnull, 'w')  

    global metabolomics_global
    metabolomics_global = m_df

def scale_metabolomics(metabolomics,fva_dir='fva/'):
    print('Scaling metabolomics...')

    raw = _load_dataframe(metabolomics)

    if not os.path.exists(fva_dir):
        raise Exception('Could not find fva directory at %s'%fva_dir)

    fva_dfs = [_load_dataframe(fva_dir+f) for f in os.listdir(fva_dir)]
    if len(fva_dfs)==0:
        raise Exception('No FVA results found')

    print('Found FVA results for %s out of %s total samples!'%(len(fva_dfs),len(raw.columns)))

    mins = pd.concat([df['min'] for df in fva_dfs],axis=1)
    maxs = pd.concat([df['max'] for df in fva_dfs],axis=1)

    min_max = pd.concat([mins.min(axis=1,numeric_only=True).rename('min'),maxs.max(axis=1,numeric_only=True).rename('max')],axis=1).T.to_dict()

    if len(fva_dfs)/len(raw.columns) < 0.5:
        print(f'!!! Using FVA results for less than 50% of all samples, results might not be as good as they could be!\n')

    print('Scaling %s metabolites...'%len(raw.index))

    missing_metabs = [m for m in list(raw.index) if m not in min_max]
    if len(missing_metabs)>0:
        print('Skipping %s metabolites due to them being missing from models, these metabolites will not be included in formatted metabolomics...\n'%len(missing_metabs)) 
        raw.drop(missing_metabs,inplace=True)

    def _scale_row(x):
        metab = x.name
        log_row = np.log10(x)

        metab_min = min_max[metab]['min']
        metab_max = min_max[metab]['max']
        scaled = ((metab_max-metab_min)*(log_row-log_row.min())/(log_row.max()-log_row.min()))+metab_min
        return scaled
    
    scaled = raw.apply(_scale_row,axis=1)
    scaled.sort_index(inplace=True)

    if len(scaled.index) == 0:
        raise Exception('Scaled metabolomics is empty. Please check FVA files to make sure everything is correct!')

    return scaled

def _add_correlation_objective(model, flux_map):
    from optlang.interface import Objective

    obj_expression = None

    flux_map = {k:v for k,v in flux_map.items() if k in model.variables}
    for f_id, flux in flux_map.items():
        forward_var = model.variables[f_id]
        reverse_var = model.variables[_get_reverse_id(f_id)]
        net = forward_var-reverse_var

        squared_diff=(net-flux)**2
        obj_expression = squared_diff if obj_expression is None else obj_expression + squared_diff
    try:
        model.objective = Objective(obj_expression,direction="min")
        model.update()
    except Exception as e:
        raise Exception('Failed to add mseFBA objective to model- %s'%e)

def _load_dataframe(m):
    if isinstance(m,str):
        if not os.path.exists(m):
            raise Exception('Tried to load dataframe from path that does not exist- %s'%m)
        return pd.read_csv(m,index_col=0)
    elif isinstance(m, pd.DataFrame):
        return m
    else:
        raise Exception('_load_dataframe can only take a string or dataframe, received %s'%type(m))

def combine_mse_solutions(mse_dir='mse/',out_file='mse_sol.csv',write_to_file=True):
    dfs = [_load_dataframe(mse_dir+f) for f in os.listdir(mse_dir)]
    combined = pd.concat(dfs,axis=1)
    combined.sort_index(axis=1,inplace=True)
    combined.sort_index(axis=0,inplace=True)
    if write_to_file:
        combined.to_csv(out_file)
    return combined

def evaluate_mse_results(metabolomics,comparison_df=None,mse_dir='mse/',map_labels=True):
    if comparison_df is None:
        comparison_df = combine_mse_solutions(mse_dir=mse_dir,write_to_file=False)
    
    comparison_df = _load_dataframe(comparison_df)
    metabolomics = _load_dataframe(metabolomics)

    if map_labels:
        conversion = _load_dataframe('sample_label_conversion.csv').conversion.to_dict()
        metabolomics.rename(conversion,axis='columns',inplace=True)

    sp_r = comparison_df.corrwith(metabolomics,method='spearman',axis=1)
    pr_r = comparison_df.corrwith(metabolomics,method='pearson',axis=1)
    combined_dataframe = pd.concat([sp_r,pr_r],axis=1)
    combined_dataframe.columns=['spearman','pearson']

    print('Avg Spearman- '+str(sp_r.mean()))
    print('Avg Pearson- '+str(pr_r.mean()))

    return combined_dataframe

def compute_nmpcs(fva_dir='fva/',out_file='nmpc_sol.csv',write_to_file=True):
    fva_files = os.listdir(fva_dir)
    print('Computing NMPCs using fva results for %s samples...'%len(fva_files))
    
    nmpcs = {}
    for f in fva_files:
        s = _load_dataframe(fva_dir+f)
        label = f.split('.csv')[0]

        nmpcs[label]=(s['min']+s['max']).to_dict()

    nmpcs = pd.DataFrame(nmpcs)
    nmpcs.sort_index(axis=1,inplace=True)
    nmpcs.sort_index(axis=0,inplace=True)
    if write_to_file:
        nmpcs.to_csv(out_file)
    return nmpcs
