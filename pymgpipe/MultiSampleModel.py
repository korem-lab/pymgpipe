
import os
from .optlang_util import _get_solver_interface,_load_cplex_model,_load_gurobi_model, suppress_stdout, load_model
import optlang 
from optlang.interface import *
import tqdm

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
        self.configuration = baseclass.Configuration(
            problem=self
        )

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
