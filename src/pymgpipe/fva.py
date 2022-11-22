import os
import pandas as pd
from gurobipy import *
from multiprocessing import Pool
from functools import partial
import tqdm
import gc
import sys
from .utils import *
from optlang.interface import Objective
from pathlib import Path
import numpy as np

def regularFVA(
    model=None,
    reactions=None,
    regex=None,
    ex_only=True,
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    write_to_file=False,
    out_dir='fva/',
    parallel=True):
    gc.enable()

    model = load_model(path=model,solver=solver)

    if reactions is None and regex is None and ex_only is True:
        regex = Constants.EX_REGEX

    reactions_to_run = [r.name for r in get_reactions(model,reactions,regex)]
        
    result_df = []
    if write_to_file:
        Path(out_dir).mkdir(exist_ok=True)
        out_file=out_dir+'%s.csv'%model.name

        if os.path.exists(out_file):
            result_df = pd.read_csv(out_file)
            result_df=  result_df[~result_df['id'].isnull()][~result_df.id.duplicated(keep='first')]
            metabs_to_skip = list(result_df.id)
            print('\nFound existing file, skipping %s metabolites...'%len(metabs_to_skip))

            reactions_to_run = [r for r in reactions_to_run if r not in metabs_to_skip] 
            result_df = result_df.to_dict('records')

    if len(reactions_to_run) == 0:
        print('---Finished! No reactions left to run---')
        return
    
    if parallel is False:
        threads=1
        parallel=False
    else:
        threads = os.cpu_count() if threads == -1 or threads > os.cpu_count() else threads
        threads = min(threads,len(reactions_to_run))

        parallel = False if threads <= 1 else parallel

    split_reactions = np.array_split(reactions_to_run,threads)
    print('Starting parallel FVA with %s chunks on %s threads'%(threads,len(split_reactions)))

    _func = _optlang_worker
    if parallel:
        p = Pool(processes=threads,initializer=partial(_pool_init,model))
        res = p.imap(_func, split_reactions)
    else:
        _pool_init(model)
        res = map(_func, split_reactions)

    for result in res:
        result_df = result_df + result

        out_df = pd.DataFrame.from_records(result_df,index='id')
        out_df.sort_index(inplace=True)
        if write_to_file:
            out_df.to_csv(out_file) 
    
    try:
        p.close()
        p.join()
    except:
        pass

    return out_df

# works for both cplex and gurobi
def _optlang_worker(metabolites):
    global global_model

    result = []
    for m in metabolites:
        forward_var = global_model.variables[m]
        reverse_var = global_model.variables[get_reverse_id(m)]
        net = forward_var - reverse_var
        
        global_model.objective = Objective(net,direction='max')
        max_sol = solve_model(model=global_model,reactions=[m]).to_dict()[global_model.name][m]

        global_model.objective = Objective(net,direction='min')
        min_sol = solve_model(model=global_model,reactions=[m]).to_dict()[global_model.name][m]

        result.append({'id':m,'min':min_sol,'max':max_sol})
    return result

def _pool_init(sample_model):
    sys.stdout = open(os.devnull, 'w')  

    global global_model
    global_model = sample_model

    global global_problem
    global_problem = sample_model.problem

