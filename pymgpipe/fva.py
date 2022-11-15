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

    model = load_model(path=model,solver=solver) if isinstance(model,str) else model

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
    
    threads = os.cpu_count() if threads == -1 else threads
    threads = min(threads,len(reactions_to_run))

    parallel = False if threads <= 1 else parallel

    # _func = _gurobi_worker if solver=='gurobi' else _optlang_worker <- might speed things up?
    _func = _optlang_worker
    if parallel:
        print('Starting parallel FVA on %s with %s reactions...'%(model.name,len(reactions_to_run)))
        p = Pool(processes=threads,initializer=partial(_pool_init,model))
        res = p.imap(_func, reactions_to_run)
    else:
        print('Starting serial FVA on %s with %s reactions...'%(model.name,len(reactions_to_run)))
        _pool_init(model)
        res = map(_func, reactions_to_run)

    for result in tqdm.tqdm(res,total=len(reactions_to_run)):
        result_df.append(result)

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
def _optlang_worker(metabolite):
    global global_model

    forward_var = global_model.variables[metabolite]
    reverse_var = global_model.variables[get_reverse_id(metabolite)]
    net = forward_var - reverse_var
    
    global_model.objective = Objective(net,direction='max')
    max_sol = solve_model(model=global_model).to_dict()[global_model.name][metabolite]

    global_model.objective = Objective(net,direction='min')
    min_sol = solve_model(model=global_model).to_dict()[global_model.name][metabolite]

    return {'id':metabolite,'min':min_sol,'max':max_sol}


def _pool_init(sample_model):
    sys.stdout = open(os.devnull, 'w')  

    global global_model
    global_model = sample_model

    global global_problem
    global_problem = sample_model.problem

