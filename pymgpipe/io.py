import os
import cobra
import os.path as path
import pickle
import optlang

class UnsupportedSolverException(Exception):
    def __init__(self, msg='Unrecognized solver. Supported solvers include gurobi and cplex', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)

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
def load_cobra_model(file):
    _, ext = path.splitext(file)
    read_func = _read_funcs[ext]
    model = read_func(file)
    
    return model

# Loads either LP file or cobra file and returns optlang model representing underlying optimization problem
def load_model(path, solver='gurobi'):
    if not os.path.exists(path):
        raise Exception('Could not find model at %s'%path)

    print('Loading model from %s...'%path)
    try:
        if solver == 'gurobi':
            model = _load_gurobi_model(path)
        elif solver == 'cplex':
            model = _load_cplex_model(path)
        else:
            raise UnsupportedSolverException
        optlang_model = _model_interfaces[solver].Model(problem=model,name=path.split('/')[-1].split('.')[0])
    except:
        optlang_model = load_cobra_model(path).solver

    return optlang_model

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