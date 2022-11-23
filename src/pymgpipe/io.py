import os
import sys
import cobra
import os.path as path
import pickle
import optlang
import logging
from contextlib import contextmanager

class UnsupportedSolverException(Exception):
    def __init__(self, msg='Unrecognized solver. Supported solvers include gurobi and cplex', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

_read_funcs = {
    ".xml": cobra.io.read_sbml_model,
    ".gz": cobra.io.read_sbml_model,
    ".mat": cobra.io.load_matlab_model,
    ".json": cobra.io.load_json_model,
    ".pickle": lambda fn: pickle.load(open(fn, "rb")),
}

_model_interfaces = {
    "cplex": optlang.cplex_interface,
    "gurobi": optlang.gurobi_interface
}

# Loads cobra file and returns cobrapy model
# RETURNS- cobra model
def load_cobra_model(file,solver='gurobi'):
    _, ext = path.splitext(file)
    read_func = _read_funcs[ext]
    try:
        with suppress_stdout():
            model = read_func(file)
    except:
        raise Exception('Error reading cobra model at %s'%file)
    model.name=file.split('/')[-1].split('.')[0]
    model.solver=solver
    
    return model

# Loads either LP file or cobra file and returns optlang model representing underlying optimization problem
# RETURNS- optlang model

def load_model(path, solver='gurobi'):
    if isinstance(path,cobra.Model):
        path.solver.name = path.name
        return path.solver
    elif isinstance(path,optlang.interface.Model):
        return path 
    elif not isinstance(path,str):
        raise Exception('Expected string, received %s'%type(path))

    if not os.path.isfile(path):
        raise Exception('Could not find model at %s'%path)

    print('Loading model from %s...'%path)
    try:
        if solver == 'gurobi':
            model = _load_gurobi_model(path)
        elif solver == 'cplex':
            model = _load_cplex_model(path)
        else:
            raise UnsupportedSolverException
        try:
            optlang_model = _model_interfaces[solver].Model(problem=model,name=path.split('/')[-1].split('.')[0])
        except:
            raise Exception('Unable to create optlang %s model, try switching solver'%solver.upper())
    except:
        optlang_model = load_cobra_model(path,solver)
        optlang_model.solver.name=path.split('/')[-1].split('.')[0]
        return optlang_model.solver


    return optlang_model

def write_lp_problem(model,out_file=None,compress=True,force=True):
    out_file = './'+model.name+'.xml' if out_file is None else out_file
    if compress:
        out_file = out_file + '.gz'
    if not force and os.path.basename(out_file).split('.')[0] in [f.split('.')[0] for f in os.listdir(os.path.dirname(out_file))]:
        print('Model already exists!')
        return

    # some computers cant compress to gz
    try:
        model.solver.problem.write(out_file)
    except:
        if compress:
            logging.warn('Could not write LP to compressed file format, trying .7z extension')
            model.solver.problem.write(out_file.replace('.gz','.7z'))

def _load_cplex_model(path):
    try:
        import cplex
        with suppress_stdout():
            return cplex.Cplex(path)
    except Exception as e:
        raise Exception('Provided model is not a valid CPLEX model')

def _load_gurobi_model(path):
    try:
        import gurobipy
        with suppress_stdout():
            return gurobipy.read(path)
    except Exception as e:
        raise Exception('Provided model is not a valid GUROBI model- %s'%path)  