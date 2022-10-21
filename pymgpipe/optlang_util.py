import os
import optlang
import re
import pandas as pd
import sys
from contextlib import contextmanager
import warnings
import logging

warnings.filterwarnings("ignore")

class Constants:
    EX_REGEX = '^EX_.*_m$'
    EX_REGEX_MULTI_SAMPLE = '^EX_.*_m_.*$'

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

class UnsupportedSolverException(Exception):
    def __init__(self, msg='Unrecognized solver. Supported solvers include gurobi and cplex', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)

class InfeasibleModelException(Exception):
    pass

def load_model(path, solver='gurobi'):
    if not os.path.exists(path):
        raise Exception('Could not find model at %s'%path)

    if solver == 'gurobi':
        model = _load_gurobi_model(path)
    elif solver == 'cplex':
        model = _load_cplex_model(path)
    else:
        raise UnsupportedSolverException
    interface = _get_solver_interface(solver)
    optlang_model = interface.Model(problem=model,name=path.split('/')[-1].split('.')[0])
    return optlang_model

def solve_model(
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

    if regex is None and reactions is None and ex_only:
        regex = Constants.EX_REGEX
    fluxes = _get_fluxes_from_model(model,threshold=flux_threshold,regex=regex,reactions=reactions)
    res = pd.DataFrame({model.name:fluxes})
    del model
    return res

def _get_fluxes_from_model(model,reactions=None,regex=None,threshold=1e-5):
    fluxes = {}

    for forward in get_reactions(model,reactions,regex):
        r_id = get_reverse_id(forward.name)
        if r_id not in model.variables:
            continue

        reverse = model.variables[r_id]

        flux = float(forward.primal-reverse.primal)
        flux = 0 if flux == -0.0 else flux
        flux = flux if abs(flux)>threshold else 0
        fluxes[forward.name]=flux
    return fluxes

def get_reactions(model,reactions=None,regex=None):
    r = []
    if reactions is not None and len(reactions)>0:
        if isinstance(reactions[0],optlang.gurobi_interface.Variable) or isinstance(reactions[0],optlang.cplex_interface.Variable):
            r = [k for k in reactions if k.name in model.variables]
        elif isinstance(reactions[0],str):
            r = [model.variables[k] for k in reactions if k in model.variables]
        else:
            raise Exception('List of reactions need to be either IDs or reaction variables! Received- %s'%type(reactions[0]))
    elif regex is not None:
        try:
            re.compile(regex)
            r = [k for k in model.variables if re.match(regex,k.name) and 'reverse' not in k.name]
        except re.error:
            raise Exception('Invalid regex- %s'%regex)
    else:
        r = [k for k in model.variables if 'reverse' not in k.name]
    if len(r) == 0:
        logging.warn('Returning 0 reactions from model!')
    return r

def constrain_reactions(model, flux_map, threshold=0.0):
    if isinstance(flux_map, pd.Series):
        flux_map = flux_map.to_dict()
    flux_map = {k:v for k,v in flux_map.items() if k in model.variables}
    for f_id, flux in flux_map.items():
        forward_var = model.variables[f_id]
        reverse_var = model.variables[get_reverse_id(f_id)]

        if flux > 0:
            forward_var.set_bounds(flux-threshold,flux+threshold)
            reverse_var.set_bounds(0,0)
        if flux < 0:
            reverse_var.set_bounds(-flux-threshold,-flux+threshold)
            forward_var.set_bounds(0,0)
        elif flux == 0:
            forward_var.set_bounds(0,threshold)
            reverse_var.set_bounds(0,threshold)
    model.update()
    return list(flux_map.keys())
 
def set_objective(model, obj_expression, direction='min'):
    try:
        model.objective = optlang.Objective(obj_expression,direction=direction)
        model.update()
        logging.info('Set model objective!')
    except Exception as e:
        raise Exception('Failed to add objective to %s- %s\n%s'%(model.name,e,obj_expression))

def get_reverse_id(id):
    import hashlib
    if not isinstance(id, str):
        try:
            id = id.name
        except:
            raise Exception('get_reverse_id must take either string ID or optlang.Variable')

    if re.match('.*_mc.*',id) is not None:
        id, sample_num = id.split('_mc')
        sample_id = '_mc'+sample_num
        return "_".join(
            (id, "reverse", hashlib.md5(id.encode("utf-8")).hexdigest()[0:5])
        ) + sample_id
    else:
        return  "_".join((id, "reverse", hashlib.md5(id.encode("utf-8")).hexdigest()[0:5]))

def get_reverse_var(model, v):
    return model.variables[get_reverse_id(v)]

def _load_cplex_model(path):
    try:
        import cplex
        return cplex.Cplex(path)
    except Exception as e:
        raise Exception('Provided model is not a valid CPLEX model')

def _load_gurobi_model(path):
    try:
        import gurobipy
        return gurobipy.read(path)
    except Exception as e:
        raise Exception('Provided model is not a valid GUROBI model- %s'%path)  

def _get_solver_interface(str):
    if str == 'gurobi':
        return optlang.gurobi_interface
    elif str == 'cplex':
        return optlang.cplex_interface
    raise UnsupportedSolverException

def _run_cplex_checklist():
    from docplex.mp.check_list import run_docplex_check_list
    run_docplex_check_list()
