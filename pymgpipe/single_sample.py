import os
from .optlang_util import _get_solver_interface,_load_cplex_model,_load_gurobi_model, suppress_stdout, load_model, get_reactions, InfeasibleModelException, Constants, get_reverse_id
import optlang 
from optlang.interface import *
import tqdm

from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np
import re

def _l_model(m_path,solver='gurobi'):
    m = load_model(m_path,solver=solver)
    return (m,m_path.split('/')[-1].split('.')[0])

class SingleSampleModel(object):
    def __init__(self, name="", solver='gurobi', problem=None, verbosity=0, presolve=True, method='auto'):
        self.name = name
        self.solver = solver

        baseclass = _get_solver_interface(solver)
        self.__class__ = type(self.__class__.__name__,
                            (baseclass.Model, object),
                            dict(self.__class__.__dict__))

        if problem is not None:
            print('Loading from existing file...')
            name = problem.split('/')[-1].split('.')[0]
            if solver == 'gurobi':
                problem = _load_gurobi_model(problem)
            elif solver == 'cplex':
                problem = _load_cplex_model(problem)
            super(self.__class__, self).__init__(problem=problem,name=name)
        else:
            super(self.__class__, self).__init__()
        self.configuration = baseclass.Configuration(
            problem=self,
            lp_method=method,
            qp_method=method,
            presolve=presolve,
            verbosity=verbosity
        )

    def solve(
        self,
        regex=None,
        reactions=None,
        flux_threshold=1e-5,
        ex_only=True):
    
        try:
            self.optimize()
        except Exception as e:
            raise Exception(f"Error when optimizing model, make sure model is valid- ${e}")

        if self.status == 'infeasible':
            raise InfeasibleModelException('%s is infeasible!'%self.name)

        if regex is None and reactions is None and ex_only:
            regex = Constants.EX_REGEX
        fluxes = self._get_fluxes(threshold=flux_threshold,regex=regex,reactions=reactions)
        res = pd.DataFrame({self.name:fluxes})
        return res

    def _get_fluxes(self,reactions=None,regex=None,threshold=1e-5):
        fluxes = {}

        for forward in get_reactions(self,reactions,regex):
            r_id = get_reverse_id(forward.name)
            if r_id not in self.variables:
                continue

            reverse = self.variables[r_id]

            flux = float(forward.primal-reverse.primal)
            flux = 0 if flux == -0.0 else flux
            flux = flux if abs(flux)>threshold else 0
            fluxes[forward.name]=flux
        return fluxes

    def get_reactions(self,reactions=None,regex=None):
        return get_reactions(self,reactions=reactions,regex=regex)

    def write_lp(self,out=None):
        if out is None:
            out = str(hash(self)) + '.lp'
            print('Writing to %s'%out)
        self.problem.write(out)
