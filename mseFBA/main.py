import os

import pandas as pd
from gurobipy import *

from multiprocessing import Pool
from functools import partial
import tqdm
import gc

from .optlang_util import InfeasibleModelException, load_model, solve_model, _add_correlation_objective, _get_exchange_reactions
import numpy as np
import sys

def mseFBA(
    metabolomics_file,
    files=[],
    file_dir='problems/',
    ex_only=True,
    zero_unmapped_metabolites=False,
    threads=int(os.cpu_count()/2),
    out_file=None,
    solver='gurobi',
    verbosity=0,
    presolve=True,
    threshold=1e-5,
    scale=True
):
    gc.enable()

    try:
        model_files = files if files else [file_dir+m for m in os.listdir(file_dir)]
    except:
        raise Exception('Please pass in a valid model directory using the \'dir\' parameter or an explicit set of files using the \'files\' parameter')

    if not os.path.exists(metabolomics_file):
        raise Exception('Metabolomics files does not exist.')

    result = pd.DataFrame()
    threads = os.cpu_count() if threads == -1 else min(threads,len(model_files))

    print('Solver- %s'%solver.upper())
    print('Presolve- %s'%str(presolve).upper())
    print('Verbosity- %s'%verbosity)
    print('Zero unmapped metabolites- %s'%str(zero_unmapped_metabolites).upper())
    print('Scale metabolomics- %s'%str(scale).upper())
    print('------------------------------------------')
    print('Running mseFBA on %s samples using %s threads...'%(len(model_files),threads))
    
    if scale:
        metabolomics_df = scale_metabolomics(metabolomics_file)
    else:
        metabolomics_df = pd.read_csv(metabolomics_file,index_col=0)

    p = Pool(processes=threads,initializer=partial(_pool_init,metabolomics_df))
    _func = partial(
        _mseFBA_worker,
        ex_only,
        zero_unmapped_metabolites,
        solver,
        verbosity,
        presolve,
        threshold
    )
    infeasible_models = []
    for sample_res in tqdm.tqdm(p.imap(_func, model_files),total=len(model_files)):
        gc.collect()

        if isinstance(sample_res,str):
            infeasible_models.append(sample_res)
            continue

        result = pd.concat([result,sample_res],axis=1)
        result.sort_index(axis=1,inplace=True) # sort columns
        result.sort_index(axis=0,inplace=True) # sort index
        if out_file is not None:
            result.to_csv(out_file)
    p.close()
    p.join()

    print('Finished mseFBA!')
    if len(infeasible_models) > 0:
        print('!!! Some models were infeasible and could not be solved-\n')
        print(infeasible_models)

    return result

def _mseFBA_worker(ex_only, zero_unmapped_metabolites, solver, verbosity, presolve, threshold, model_file):
    model = load_model(model_file,solver)
    
    metabolite_dict = metabolomics_global[model.name].to_dict()
    _add_correlation_objective(model,metabolite_dict)

    if zero_unmapped_metabolites:
        ex_reactions = _get_exchange_reactions(model)
        unmapped_map = {k:0 for k in ex_reactions if k not in metabolite_dict}
        _add_correlation_objective(model,unmapped_map)

    try:
        solution = solve_model(
            model=model,
            ex_only=ex_only,
            verbosity=verbosity,
            presolve=presolve,
            method='barrier',
            flux_threshold=threshold
        )
        del model
        return solution
    except InfeasibleModelException:
        del model
        return model_file

def _pool_init(df):
    sys.stdout = open(os.devnull, 'w')  

    global metabolomics_global
    metabolomics_global = df

def scale_metabolomics(metabolomics_path,force=False):
    print('Scaling metabolomics at path %s...'%metabolomics_path)

    if not os.path.exists(metabolomics_path):
        raise Exception('Could not find metabolomics file at specified path- %s'%metabolomics_path)

    out_path = metabolomics_path if metabolomics_path.endswith('_scaled.csv') else metabolomics_path.split('/')[-1].split('.csv')[0]+'_scaled.csv'
    if os.path.exists(out_path) and not force:
        print('Found existing scaled metabolomics, skipping this step. If you want to repeat, delete %s.'%out_path)
        return pd.read_csv(out_path,index_col=0)

    raw = pd.read_csv(metabolomics_path,index_col=0)

    if not os.path.exists('fva/'):
        raise Exception('Could not find fva/ file directory. Please place FVA sample results in folder titled \'fva\'')

    fva_dfs = [pd.read_csv('fva/%s'%f,index_col=0) for f in os.listdir('fva/')]
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
        print('!!! Skipping %s metabolites due to them being missing from models, these metabolites will not be included in formatted metabolomics...\n'%len(missing_metabs)) 
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

    scaled.to_csv(out_path)
    return scaled

