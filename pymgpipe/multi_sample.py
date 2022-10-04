from .optlang_util import *
import tqdm
import sys
from multiprocessing import Pool
from functools import partial
from .optlang_util import _get_reverse_id, _get_solver_interface, _load_gurobi_model, _load_cplex_model
import optlang
from optlang.interface import *

import os

def l_model(m_path,solver='gurobi'):
    m = load_model(m_path,solver=solver)
    return (m,'_'+m_path.split('/')[-1].split('.')[0])

def combine_samples(
    models,
    out_file,
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    parallel=False
    ):
    
    if isinstance(models,str):
        models = [models+m for m in os.listdir(models)]
    combined_model = optlang.Model()
    obj_expr = None

    threads = os.cpu_count()-1 if threads == -1 else threads
    threads = min(threads,len(models))

    print('Combining %s models using %s threads...'%(len(models),threads))
    sys.stdout = open(os.devnull, 'w')  

    _func = partial(l_model,solver=solver)
    if parallel:
        p = Pool(processes=threads)
        res = p.imap(_func, models)
    else:
        res = map(_func, models)

    for res in tqdm.tqdm(res,total=len(models)):
        curr,label=res
        for v in curr.variables:
            v.name = v.name + label
        for c in curr.constraints:
            c.name = c.name + label

        obj_expr = curr.objective.expression if obj_expr is None else obj_expr + curr.objective.expression

        combined_model._add_variables(curr.variables)
        combined_model._add_constraints(curr.constraints)
    sys.stdout = sys.__stdout__
    combined_model.objective = optlang.Objective(obj_expr)
    combined_model.update()
    combined_model.problem.write(out_file)

def solve_multi_sample_model(
    model,
    regex=None,
    reactions=None,
    solver='gurobi',
    verbosity=0,
    presolve=True,
    method='auto',
    flux_threshold=1e-5,
    ex_only=True):
    if isinstance(model,str):
        model = load_model(model, solver)
    model.configuration.verbosity=verbosity
    model.configuration.presolve=presolve
    model.configuration.lp_method=method
    model.configuration.qp_method=method

    model.update()
    
    try:
        model.optimize()
    except Exception as e:
        raise Exception(f"Error when optimizing model, make sure model is valid- ${e}")

    if model.status == 'infeasible':
        raise InfeasibleModelException('%s is infeasible!'%model.name)

    regex = '^EX_.*_m_.*$' if ex_only else regex
    fluxes = _get_fluxes_from_model(model,threshold=flux_threshold,regex=regex,reactions=reactions)
    return pd.DataFrame(fluxes)

def _get_fluxes_from_model(model,reactions=None,regex=None,threshold=1e-5):
    fluxes = {}

    for forward in get_reactions(model,reactions,regex):
        sample_id = 'mc'+forward.name.split('_mc')[1]
        r_id = _get_reverse_id(forward.name, multi_sample=True)
        if r_id not in model.variables:
            continue
        reverse = model.variables[r_id]

        flux = float(forward.primal-reverse.primal)
        flux = 0 if flux == -0.0 else flux
        flux = flux if abs(flux)>threshold else 0
        if sample_id not in fluxes:
            fluxes[sample_id] = {}
        fluxes[sample_id][forward.name.split('_mc')[0]]=flux
    return fluxes

class MultiSampleModel(object):
    def __init__(self, name="", solver='gurobi', problem=None):
        self.name = name
        self.samples = []
        self.solver = solver

        baseclass = _get_solver_interface(solver)
        self.__class__ = type(self.__class__.__name__,
                            (baseclass.Model, object),
                            dict(self.__class__.__dict__))

        if problem is not None:
            print('Loading from existing file...')
            if solver == 'gurobi':
                problem = _load_gurobi_model(problem)
            elif solver == 'cplex':
                problem = _load_cplex_model(problem)
            super(self.__class__, self).__init__(problem=problem)
            self.samples = list(set(m.name.split('_')[-1] for m in self.variables))
        else:
            super(self.__class__, self).__init__()
        self.configuration = baseclass.Configuration()

    def add_sample(self, model, name=None):
        if isinstance(model,str):
            with suppress_stdout():
                model = load_model(model)
        name = model.name if name is None else name
        if name in self.samples:
            print('Sample %s already exists in model!'%name)
            return
       
        print('Adding %s to model...'%name)

        label = '_'+name
        for v in model.variables:
            v.name = v.name + label
        for c in model.constraints:
            c.name = c.name + label

        self._add_variables(model.variables)
        self._add_constraints(model.constraints)
    
        self.objective = optlang.Objective(self.objective.expression + model.objective.expression, direction='max')
        self.samples.append(name)
        del model

    def add_samples(self, models):
        if isinstance(models,str):
            models = [models+m for m in os.listdir(models)]
        print('Adding %s samples to model...'%len(models))
        for m in tqdm.tqdm(models):
            with suppress_stdout():
                self.add_sample(m)
     
    def write_lp(self,out=None):
        if out is None:
            out = str(hash(self)) + '.lp'
            print('Writing to %s'%out)
        self.problem.write(out)
