import os
import sys
import cobra
import os.path as path
import pickle
import optlang
import logging
from contextlib import contextmanager
from .sbml import write_sbml_model


class UnsupportedSolverException(Exception):
    def __init__(
        self,
        msg="Unrecognized solver. Supported solvers include gurobi and cplex",
        *args,
        **kwargs
    ):
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

_write_funcs = {
    ".xml": write_sbml_model,
    ".gz": write_sbml_model,
    ".mat": cobra.io.save_matlab_model,
    ".json": cobra.io.save_json_model,
    ".pickle": lambda fn: pickle.dump(open(fn, "wb")),
}


def show_availible_solvers():
    return [k.lower() for k, v in optlang.available_solvers.items() if v]


def _get_optlang_interface(solver):
    available_solvers = show_availible_solvers()
    if solver not in available_solvers:
        raise Exception(
            "Provided solver %s not available. Available solvers are %s"
            % (solver, available_solvers)
        )
    elif solver == "gurobi":
        return optlang.gurobi_interface
    elif solver == "cplex":
        return optlang.cplex_interface
    else:
        raise Exception(
            "Provided solver %s is unsupported. Solver must either be `gurobi` or `cplex`."
            % solver
        )


# Loads cobra file and returns cobrapy model
# RETURNS- cobra model
def load_cobra_model(file, solver="gurobi"):
    _, ext = path.splitext(file)
    read_func = _read_funcs[ext]
    try:
        with suppress_stdout():
            model = read_func(file)
    except Exception:
        raise Exception("Error reading cobra model at %s" % file)
    model.name = file.split("/")[-1].split(".")[0]
    model.solver = solver

    return model


def write_cobra_model(model, file):
    _, ext = path.splitext(file)
    write_func = _write_funcs[ext]
    try:
        with suppress_stdout():
            write_func(model, file)
    except Exception:
        raise Exception("Error writing cobra model to %s" % file)


# Loads either LP file or cobra file and returns optlang model representing underlying optimization problem
# RETURNS- optlang model
def load_model(path, solver="gurobi"):
    if isinstance(path, cobra.Model):
        path.solver.name = path.name
        return path.solver
    elif isinstance(path, optlang.interface.Model):
        return path
    elif not isinstance(path, str):
        raise Exception("Expected string, received %s" % type(path))

    if not os.path.isfile(path):
        raise Exception("Could not find model at %s" % path)

    print("Loading model from %s..." % path)
    try:
        if solver == "gurobi":
            model = _load_gurobi_model(path)
        elif solver == "cplex":
            model = _load_cplex_model(path)
        else:
            raise UnsupportedSolverException
        try:
            optlang_model = _get_optlang_interface(solver).Model(
                problem=model, name=path.split("/")[-1].split(".")[0]
            )
        except Exception:
            raise Exception(
                "Unable to create optlang %s model, try switching solver"
                % solver.upper()
            )
    except Exception:
        optlang_model = load_cobra_model(path, solver)
        optlang_model.solver.name = path.split("/")[-1].split(".")[0]
        return optlang_model.solver

    return optlang_model


def write_lp_problem(model, out_file=None, compress=True, force=True):
    out_file = "./" + model.name + ".xml" if out_file is None else out_file
    if compress:
        out_file = out_file + ".gz"
    if not force and os.path.basename(out_file).split(".")[0] in [
        f.split(".")[0] for f in os.listdir(os.path.dirname(out_file))
    ]:
        print("Model already exists!")
        return

    # some computers cant compress to gz
    try:
        model.solver.problem.write(out_file)
    except Exception:
        if compress:
            logging.warn(
                "Could not write LP to compressed file format, trying .7z extension"
            )
            model.solver.problem.write(out_file.replace(".gz", ".7z"))


def _load_cplex_model(path):
    try:
        import cplex

        with suppress_stdout():
            return cplex.Cplex(path)
    except Exception:
        raise Exception("Provided model is not a valid GUROBI model- %s" % path)


def _load_gurobi_model(path):
    try:
        import gurobipy

        with suppress_stdout():
            return gurobipy.read(path)
    except Exception:
        raise Exception("Provided model is not a valid GUROBI model- %s" % path)
