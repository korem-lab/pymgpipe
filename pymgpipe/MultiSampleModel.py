
import os
from .optlang_util import _get_solver_interface,_load_cplex_model,_load_gurobi_model, suppress_stdout, load_model, get_reactions, InfeasibleModelException, Constants
import optlang 
from optlang.interface import *
import tqdm

from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np

def _l_model(m_path,solver='gurobi'):
    m = load_model(m_path,solver=solver)
    return (m,m_path.split('/')[-1].split('.')[0])

class MultiSampleModel(object):
    def __init__(self, name="", solver='gurobi', problem=None, samples=None, verbosity=0, presolve=True, method='auto'):
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
        self.configuration = baseclass.Configuration(
            problem=self,
            lp_method=method,
            qp_method=method,
            presolve=presolve,
            verbosity=verbosity
        )
        if samples is not None:
            self.add_samples(samples,threads=-1)

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

    def add_samples(
        self,
        models,
        parallel=True,
        threads=int(os.cpu_count()/2)):
        if isinstance(models,str):
            if os.path.isdir(models):
                models = [models+m for m in os.listdir(models)]
            elif os.path.exists(models):
                models = [models]
            else:
                raise Exception('If passing in models parameter as a string, it must be valid model path or directory')
        threads = os.cpu_count()-1 if threads == -1 else threads
        threads = min(threads,len(models))

        parallel = False if threads == 1 else parallel
        if parallel:
            print('Adding %s samples using %s threads...'%(len(models),threads))

            _func = partial(_l_model,solver=self.solver)
            p = Pool(processes=threads)

            sys.stdout = open(os.devnull, 'w')  
            for m,name in tqdm.tqdm(p.imap(_func,models),total=len(models)):
                    self.add_sample(m,name=name)
            sys.stdout = sys.__stdout__

            p.close()
            p.join()
        else:
            print('Adding %s samples to model in series...'%len(models))

            for m in tqdm.tqdm(models):
                with suppress_stdout():
                    self.add_sample(m)
            
    def solve(
        self,
        regex=None,
        reactions=None,
        flux_threshold=1e-5,
        ex_only=True):
        self.update()

        try:
            self.optimize()
        except Exception as e:
            raise Exception(f"Error when optimizing model, make sure model is valid- ${e}")

        if self.status == 'infeasible':
            raise InfeasibleModelException('%s is infeasible!'%self.name)

        regex = Constants.EX_REGEX_MULTI_SAMPLE if ex_only else regex
        fluxes = self._get_fluxes(threshold=flux_threshold,regex=regex,reactions=reactions)
        return pd.DataFrame(fluxes)

    def compute_nmpcs(self,reactions=None,ex_only=True):
        comb = pd.DataFrame()
        if reactions is None and ex_only is True:
            reactions = list(set([m.name.split('_mc')[0] for m in get_reactions(self,regex=Constants.EX_REGEX_MULTI_SAMPLE)]))

        print('Performing FVA on %s reactions...'%len(reactions))
        for metab in tqdm.tqdm(reactions):
            metab_rxns = get_reactions(self,regex='%s_mc.*$'%metab)
            net_rxns = [f-self.variables[MultiSampleModel.get_reverse_id(f.name)] for f in metab_rxns]
            s_obj = np.sum(net_rxns)
            
            self.objective = optlang.Objective(s_obj,direction='min')
            min_sol = self.solve(reactions=metab_rxns)
            
            self.objective = optlang.Objective(s_obj,direction='max')
            max_sol = self.solve(reactions=metab_rxns)

            comb = pd.concat([comb,min_sol.add(max_sol)],axis=0)
        return comb

    def _get_fluxes(self,reactions=None,regex=None,threshold=1e-5):
        fluxes = {}
        rxns = get_reactions(self,reactions=reactions,regex=regex)
        for forward in rxns:
            sample_id = 'mc'+forward.name.split('_mc')[1]
            r_id = MultiSampleModel.get_reverse_id(forward.name)
            if r_id not in self.variables:
                continue
            reverse = self.variables[r_id]

            flux = float(forward.primal-reverse.primal)
            flux = 0 if flux == -0.0 else flux
            flux = flux if abs(flux)>threshold else 0
            if sample_id not in fluxes:
                fluxes[sample_id] = {}
            fluxes[sample_id][forward.name.split('_mc')[0]]=flux
        return fluxes

    @staticmethod
    def get_reverse_id(id):
        import hashlib
        id, sample_num = id.split('_mc')
        sample_id = '_mc'+sample_num
        return "_".join(
            (id, "reverse", hashlib.md5(id.encode("utf-8")).hexdigest()[0:5])
        ) + sample_id


    def write_lp(self,out=None):
        if out is None:
            out = str(hash(self)) + '.lp'
            print('Writing to %s'%out)
        self.problem.write(out)
