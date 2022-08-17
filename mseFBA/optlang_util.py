import os
from queue import Empty
import optlang
import gurobipy
import cplex

import re
import pandas as pd

class UnsupportedSolverException(Exception):
    def __init__(self, msg='Unrecognized solver. Supported solvers include gurobi and cplex', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)

class InfeasibleModelException(Exception):
    pass

class Constants:
    EX_REGEX = '^EX_.*_m$'

def load_model(path, solver='gurobi'):
    if not os.path.exists(path):
        raise Exception('Could not find model at %s'%path)

    if path[-4:] != '.mps':
        raise Exception('Path must point to an .mps file')

    if solver == 'gurobi':
        model = _load_gurobi_model(path)
        if model is None:
            raise Exception('Provided model is not a valid GUROBI model')  
    elif solver == 'cplex':
        model = _load_cplex_model(path)
        if model is None:
            raise Exception('Provided model is not a valid CPLEX model')
    else:
        raise UnsupportedSolverException
    
    interface = _get_solver_interface(solver)
    optlang_model = interface.Model(problem=model,name=path.split('.')[0].split('/')[-1])
    return optlang_model

def solve_model(
    path=None,
    model=None,
    ex_only=False,
    solver='gurobi',
    verbosity=0,
    presolve=True,
    method='auto',
    flux_threshold=1e-5):
    if model is None:
        model = load_model(path, solver)
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

    all_fluxes = _get_fluxes_from_model(model,threshold=flux_threshold)
    if ex_only:
        all_fluxes = {k:v for k,v in all_fluxes.items() if re.match(Constants.EX_REGEX,k)}

    return pd.DataFrame({model.name:all_fluxes})

def constrain_reactions(model, flux_map, threshold=0.0):
    flux_map = {k:v for k,v in flux_map.items() if k in model.variables}
    for f_id, flux in flux_map.items():
        forward_var = model.variables[f_id]
        reverse_var = model.variables[_get_reverse_id(f_id)]

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
    
    model.objective = Objective(obj_expression,direction="min")
    model.update()

def _get_fluxes_from_model(model,specific_reactions=None,threshold=1e-5):
    fluxes = {}
    for i in range(0,len(model.variables)-1,2):
        forward = model.variables[i]
        reverse = model.variables[i+1]

        if specific_reactions is not None and forward.name not in specific_reactions:
            continue

        flux = float(forward.primal-reverse.primal)
        flux = 0 if flux == -0.0 else flux
        flux = flux if abs(flux)>threshold else 0
        fluxes[forward.name]=flux
    return fluxes

def _get_reverse_id(id):
    import hashlib
    return "_".join(
        (id, "reverse", hashlib.md5(id.encode("utf-8")).hexdigest()[0:5])
    )

def _get_exchange_reactions(model):
    return [k.name for k in model.variables if re.match(Constants.EX_REGEX,k.name)]

def _get_all_forward_reactions(model):
    return [k.name for k in model.variables if 'reverse' not in k.name]

def _load_cplex_model(path):
    try:
        return cplex.Cplex(path)
    except Exception as e:
        return None

def _load_gurobi_model(path):
    try:
        return gurobipy.read(path)
    except Exception as e:
        raise None

def _get_solver_interface(str):
    if str == 'gurobi':
        return optlang.gurobi_interface
    elif str == 'cplex':
        return optlang.cplex_interface
    raise UnsupportedSolverException

def _run_cplex_checklist():
    from docplex.mp.check_list import run_docplex_check_list
    run_docplex_check_list()
