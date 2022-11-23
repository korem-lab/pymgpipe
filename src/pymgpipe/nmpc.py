import pandas as pd
import os
from .fva import regularFVA
from .utils import *
import cobra
from collections import namedtuple
import tqdm
from pathlib import Path

def compute_nmpcs(
    samples,
    out_dir='./',
    out_file = 'nmpcs.csv',
    objective_out_file = 'community_objectives.csv',
    fluxes_out_file = 'all_fluxes.csv',
    reactions=None,
    regex=None,
    ex_only=True,
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    parallel=True,
    diet_fecal_compartments=True,
    force=False
):
    if not os.path.isdir(out_dir):
        raise Exception('Provided out directory %s is not a valid directory!'%out_dir)
    Path(out_dir).mkdir(exist_ok=True)

    out_file = out_dir+out_file
    objective_out_file = out_dir+objective_out_file
    fluxes_out_file = out_dir+fluxes_out_file

    nmpcs = pd.DataFrame() if force else load_dataframe(out_file,return_empty=True)
    all_fluxes = pd.DataFrame() if force else load_dataframe(fluxes_out_file,return_empty=True)
    obj_values = pd.DataFrame() if force else load_dataframe(objective_out_file,return_empty=True)
    obj_values['communityBiomass'] = None if obj_values.empty else obj_values[obj_values.columns[0]]

    try:
        models = [load_model(samples)]
    except:
        models = samples if isinstance(samples,list) else [os.path.dirname(samples)+'/'+m for m in os.listdir(os.path.dirname(samples))]
    models = [f for f in models if not isinstance(f,str) or os.path.basename(f).split('.')[0] not in list(nmpcs.columns)]
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
        if m.status == 'infeasible':
            logging.warn('%s model is infeasible!'%m.name)
            continue
        
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
            try:
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
            except:
                logging.warn('Cannot solve %s model!'%m.name)
                continue
        if res is None:
            return 
        res['sample_id']=m.name
        all_fluxes = pd.concat([all_fluxes,res],axis=0)
        if diet_fecal_compartments:
            metabs = [m for m in res.index.str.split('[').str[0].drop_duplicates() if not m.startswith('Diet')]
            df = {}
            for metab in metabs:
                fe = res.loc[metab+'[fe]']['max']
                d = res.loc['Diet_'+metab+'[d]']['min']
                
                df[metab.split('EX_')[1]]=d+fe
            nmpc = pd.DataFrame({m.name:df})
        else:
            nmpc = res['min']+res['max']
            nmpc.name = m.name

        nmpcs = pd.concat([nmpcs,nmpc],axis=1)
        nmpcs.to_csv(out_file)
        obj_values.to_csv(objective_out_file)
        all_fluxes.to_csv(fluxes_out_file)

    res = namedtuple('res', 'nmpc objectives fluxes')     
    return res(nmpcs,obj_values,all_fluxes)

