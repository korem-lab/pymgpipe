import os
import sys
import pandas as pd
from multiprocessing import Pool
from functools import partial
import tqdm
import gc

from pymgpipe import load_model, solve_model
from pymgpipe.optlang_util import _get_exchange_reactions, _get_reverse_id

from .metabolomics import *
from .utils import *

from optlang.interface import Objective

class paths(object):
    def __init__(self, ds_dir='./'):
        self.dataset_directory = ds_dir
        self.refresh_all()

    @property
    def dataset_directory(self):
        return self._dataset_directory

    @dataset_directory.setter
    def dataset_directory(self, value):
        self._dataset_directory = value

    @property
    def fva_directory(self):
        return self._fva_directory

    @fva_directory.setter
    def fva_directory(self, value):
        self._fva_directory = value

    @property
    def problem_directory(self):
        return self._problem_directory

    @problem_directory.setter
    def problem_directory(self, value):
        self._problem_directory = value

    @property
    def conversion_file(self):
        return self._conversion_file

    @conversion_file.setter
    def conversion_file(self, value):
        self._conversion_file = value

    def refresh_all(self):
        self.fva_directory = self.dataset_directory + 'fva/'
        self.problem_directory = self.dataset_directory + 'problems/'
        self.conversion_file = self.dataset_directory + 'sample_label_conversion.csv'

def run(
    metabolomics,
    out_file,
    dataset_dir='./',
    samples='problems/',
    fva_dir='fva/',
    conversion_file='sample_label_conversion.csv',
    ex_only=True,
    zero_unmapped_metabolites=False,
    threads=int(os.cpu_count()/2),
    solver='gurobi',
    verbosity=0,
    presolve=True,
    threshold=1e-5,
    scale=True,
    map_labels=True,
    parallelize=True,
    problem_type = '.mps',
):
    gc.enable()
    fva_dir = dataset_dir+fva_dir
    conversion_file = dataset_dir+conversion_file
 
    try:
        model_files = samples if isinstance(samples,list) else [dataset_dir+samples+m for m in os.listdir(dataset_dir+samples)]
    except:
        raise Exception('Please pass in a valid model directory or an explicit list of sample paths using the \'problems\' parameter')

    metabolomics_df = process_metabolomics(metabolomics, fva_dir, scale, map_labels, conversion_file)
    unmatched_metabolomics = [f for f in model_files if f.split('/')[-1].split(problem_type)[0] not in list(metabolomics_df.columns)]
    if len(unmatched_metabolomics) > 0:
        print('%s samples dont have associated columns in metabolomics file-\n'%len(unmatched_metabolomics))
        print([f.split('/')[-1].split(problem_type)[0] for f in unmatched_metabolomics])
        model_files = [f for f in model_files if f not in unmatched_metabolomics]
        
    solution_df = load_dataframe(out_file,return_empty=True)
    finished = list(solution_df.columns)
    if len(finished)>0:
        print('Skipping %s samples that are already finished!'%len(finished))
        model_files = [f for f in model_files if f.split('/')[-1].split(problem_type)[0] not in finished]

    if len(model_files) == 0:
        print('Finished mseFBA, no samples left to run!')
        return

    threads = os.cpu_count() if threads == -1 else threads
    threads = min(threads,len(model_files))

    print('\n---------------Parameters---------------')
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
        threshold
    )
    
    final = load_dataframe(out_file,return_empty=True)

    if parallelize:
        print('Running mseFBA on %s samples in parallel using %s threads...\n'%(len(model_files),threads))
        p = Pool(processes=threads,initializer=partial(_pool_init,metabolomics_df))
        res = p.imap(_func, model_files)
    else:
        print('Running mseFBA on %s samples in series...\n'%len(model_files))
        _pool_init(metabolomics_df)
        res = map(_func, model_files)
    
    infeasible = []
    feasible = []
    for r in tqdm.tqdm(res,total=len(model_files)):
        sample_file, solution = r
        if solution is not None:
            final = pd.concat([final,solution],axis=1)
            final.sort_index(inplace=True)
            final.to_csv(out_file)
            feasible.append(sample_file)
        else:
            infeasible.append(sample_file)
    try:
        p.close()
        p.join()
    except:
        pass

    sys.stdout = sys.__stdout__
    print('Finished mseFBA! Solved %s samples and saved to %s'%(len(feasible),out_file))
    if len(infeasible) > 0:
        print('Some models were infeasible and could not be solved-\n')
        print(list(zip(*infeasible))[0])

def _mseFBA_worker(ex_only, zero_unmapped_metabolites, solver, verbosity, presolve, threshold, model_file):
    global solution_path, metabolomics_global
    model = load_model(model_file,solver)
    if model.name not in metabolomics_global.columns:
        raise Exception('Could not find %s in metabolomics- %s'%(model_file,list(metabolomics_global.columns)))

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
    except:
        pass

    del model
    gc.collect()

    return (model_file, solution)
 
def _add_correlation_objective(model, flux_map):
    obj_expression = None

    flux_map = {k+'_'+model.name:v for k,v in flux_map.items() if k in model.variables or k+'_'+model.name in model.variables}
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

def _pool_init(m_df):
    sys.stdout = open(os.devnull, 'w')  

    global metabolomics_global
    metabolomics_global = m_df
