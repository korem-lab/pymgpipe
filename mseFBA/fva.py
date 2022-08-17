import os

import pandas as pd
from gurobipy import *

from multiprocessing import Pool
from functools import partial
import tqdm
import gc
import sys

from .optlang_util import load_model, solve_model, _get_exchange_reactions, _get_all_forward_reactions, _get_reverse_id
from optlang.interface import Objective
from pathlib import Path

def regularFVA(
    model_path, 
    ex_only = True, 
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    write_to_file=True):
    gc.enable()

    if not os.path.exists(model_path):
        raise Exception('Could not find sample at path %s'%model_path)
    
    if write_to_file:
        Path('fva').mkdir(exist_ok=True)

    model = load_model(path=model_path,solver=solver)
    reactions_to_run = _get_exchange_reactions(model) if ex_only else _get_all_forward_reactions(model)

    print('Starting FVA on %s with %s reactions...'%(model.name,len(reactions_to_run)))

    out_file='fva/%s.csv'%model.name
    result_df = []
    if os.path.exists(out_file):
        result_df = pd.read_csv(out_file)
        metabs_to_skip = result_df['id'].unique()
        print('\nFound existing file, skipping %s metabolites...'%len(metabs_to_skip))

        reactions_to_run = [r for r in reactions_to_run if r not in metabs_to_skip] 
        result_df = result_df.to_dict('records')
    
    threads = os.cpu_count() if threads == -1 else min(threads,len(reactions_to_run))

    print('Solving for %s remaining metabolites using %s threads...'%(len(reactions_to_run),threads))

    p = Pool(processes=threads,initializer=partial(_pool_init,model))
    _func = _single_metabolite_worker

    failed_metabolites = []
    out_df = pd.DataFrame()

    for result in tqdm.tqdm(p.imap(_func, reactions_to_run),total=len(reactions_to_run)):
        if isinstance(result, str):
            failed_metabolites.append(result)
        result_df.append(result)

        out_df = pd.DataFrame.from_records(result_df,index='id')
        out_df.sort_index(inplace=True)
        if write_to_file:
            out_df.to_csv(out_file) 

    p.close()
    p.join()

    return out_df

def _single_metabolite_worker(metabolite):
    global global_model

    if metabolite not in global_model.variables:
        return metabolite
    forward_var = global_model.variables[metabolite]
    reverse_var = global_model.variables[_get_reverse_id(metabolite)]
    
    global_model.objective = Objective(forward_var-reverse_var,direction='max')
    max_sol = solve_model(model=global_model).to_dict()[global_model.name][metabolite]

    global_model.objective = Objective(forward_var-reverse_var,direction='min')
    min_sol =solve_model(model=global_model).to_dict()[global_model.name][metabolite]

    return {'id':metabolite,'min':min_sol,'max':max_sol}

def _pool_init(sample_model):
    sys.stdout = open(os.devnull, 'w')  

    global global_model
    global_model = sample_model