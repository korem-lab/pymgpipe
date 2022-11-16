import pandas as pd
import os
from .fva import regularFVA
from .utils import *
import cobra
from collections import namedtuple
import tqdm


def compute_nmpcs(
    samples=[],
    out_file = 'nmpcs.csv',
    reactions=None,
    regex=None,
    ex_only=True,
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    parallel=True,
    diet_fecal_compartments=True,
    force=False
):
    objective_out_file = '/'.join(out_file.split('/')[:-1])+'/community_objectives.csv' if '/' in out_file else 'community_objectives.csv'
    
    nmpcs = pd.DataFrame() if force else load_dataframe(out_file,return_empty=True)
    obj_values = pd.DataFrame(columns=['communityObj']) if force else load_dataframe(objective_out_file,return_empty=True)

    if isinstance(samples,str) and os.path.isfile(samples):
        models = [samples]
    else:
        models = samples if isinstance(samples,list) else [samples.split('/')[0]+'/'+m for m in os.listdir(samples.split('/')[0]+'/')]
    models = [f for f in models if not isinstance(f,str) or f.split('/')[-1].split('.')[0] not in list(nmpcs.columns)]

    threads = os.cpu_count()-1 if threads == -1 else threads 

    print('Computing NMPCs on %s models...'%len(models))
    print('\n----------------Parameters----------------')
    print('Parallel- %s'%str(parallel).upper())
    print('Threads- %s'%str(threads).upper())
    print('Solver- %s'%solver.upper())
    print('Diet/fecal compartments- %s'%str(diet_fecal_compartments).upper())
    print('Exchanges only- %s'%str(ex_only).upper())
    print('Out file- %s'%str(out_file).upper())
    print('Force ovewrite- %s'%str(force).upper())
    print('------------------------------------------')

    for s in tqdm.tqdm(models,total=len(models)):
        with suppress_stdout():
            m = load_model(path=s,solver=solver)
            if not isinstance(m,optlang.interface.Model):
                raise Exception('Expected optlang.Model, received %s'%type(m))
            if not force and m.name in list(nmpcs.columns):
                continue
    
        # Solve for objective first
        if 'communityBiomass' not in m.variables:
            raise Exception('Could not find communityBiomass variable in model!')
        m.variables['communityBiomass'].set_bounds(0.4,1)
        set_objective(m,m.variables['communityBiomass'],direction='max')
        m.optimize()
        
        obj_val = round(m.objective.value,5)
        obj_values.loc[m.name]=obj_val
        if 'ObjectiveConstraint' in m.constraints:
            m.remove(m.constraints['ObjectiveConstraint'])
            m.update()
        obj_const = m.interface.Constraint(expression=m.objective.expression,lb=obj_val,ub=obj_val,name='ObjectiveConstraint')
        m.add(obj_const)
        m.update()

        # Now perform FVA under constrained objective value
        with suppress_stdout():
            res = regularFVA(
                m,
                reactions=reactions,
                regex=regex,
                ex_only=ex_only,
                solver=solver,
                threads=threads,
                parallel=parallel,
                write_to_file=False
            )
        if res is None:
            return 
        if diet_fecal_compartments:
            metabs = [m for m in res.index.str.split('[').str[0].drop_duplicates() if not m.startswith('Diet')]
            df = {}
            for metab in metabs:
                fe = res.loc[metab+'[fe]']['max']
                d = res.loc['Diet_'+metab+'[d]']['min']
                
                df[metab]=d+fe
            nmpc = pd.DataFrame({m.name:df})
        else:
            nmpc = res['min']+res['max']
            nmpc.name = m.name

        nmpcs = pd.concat([nmpcs,nmpc],axis=1)
        nmpcs.to_csv(out_file)
        obj_values.to_csv(objective_out_file)

    res = namedtuple('res', 'nmpc objectives')     
    return res(nmpcs,obj_values)

