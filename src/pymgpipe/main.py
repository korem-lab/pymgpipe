import sys
sys.path.insert(1, '../')

import pandas as pd
import os
import cobra
from .io import load_cobra_model, write_lp_problem, write_cobra_model
from pathlib import Path
from pkg_resources import resource_listdir,resource_filename

from multiprocessing import Pool
from functools import partial
import gc
import tqdm
import pickle
from .build import build
import re
from .diet import *

cobra_config = cobra.Configuration()
cobra_config.lower_bound=-1000
cobra_config.upper_bound=1000

def build_models(
    coverage_file,
    taxa_dir,
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    samples=None,
    parallel=True,
    lp_type='.mps',
    cobra_type='.xml',
    out_dir='./',
    coupling_constraints=True,
    diet_fecal_compartments=False,
    diet=None,
    compress=True,
    write_lp=True
):   
    cobra_config.solver = solver

    gc.disable()
    Path(out_dir).mkdir(exist_ok=True)

    taxa_dir = taxa_dir+'/' if taxa_dir[-1] != '/' else taxa_dir
    out_dir = out_dir+'/' if out_dir[-1] != '/' else out_dir

    model_dir=out_dir+'models/'
    problem_dir=out_dir+'problems/'
    Path(model_dir).mkdir(exist_ok=True)
    if write_lp:
        Path(problem_dir).mkdir(exist_ok=True)

    formatted = _format_coverage_file(coverage_file,taxa_dir,out_dir)
    samples_to_run = formatted.sample_id.unique()
    taxa = formatted.strain.unique()

    print('Found coverage file with %s samples and %s unique taxa'%(len(samples_to_run),len(taxa)))
   
    if samples is not None:
        samples_to_run = samples if isinstance(samples,list) else [samples]
    
    if len(samples_to_run)<=1:
        parallel=False

    threads = os.cpu_count()-1 if threads == -1 else threads
    threads = min(threads,len(samples_to_run))

    print('Building %s models...'%len(samples_to_run))
    print('\n----------------Parameters----------------')
    print('Diet/fecal compartments- %s'%str(diet_fecal_compartments).upper())
    print('Coupling constraints- %s'%str(coupling_constraints).upper())
    print('Parallel- %s'%str(parallel).upper())
    print('Threads- %s'%str(threads).upper())
    print('Solver- %s'%solver.upper())
    print('LP type- %s'%str(lp_type.split('.')[1]).upper())
    print('COBRA type- %s'%str(cobra_type.split('.')[1]).upper())
    print('compress- %s'%str(compress).upper())
    print('Output directory- %s'%str(out_dir).upper())

    _func = partial(
        _build_single_model,
        formatted,
        solver,
        model_dir,
        problem_dir,
        lp_type,
        cobra_type,
        coupling_constraints,
        diet_fecal_compartments,
        diet,
        compress,
        write_lp
    )
    
    if parallel:        
        p = Pool(threads, initializer=_mute)
        p.daemon = False

        built = list(tqdm.tqdm(p.imap(_func, samples_to_run),total=len(samples_to_run)))

        p.close()
        p.join()
    else:
        built = tqdm.tqdm(list(map(_func,samples_to_run)),total=len(samples_to_run))
    print('-------------------------------------------------------')
    print('Finished building %s models and associated LP problems!'%len(built))    

def _build_single_model(coverage_df,solver,model_dir,problem_dir,lp_type,cobra_type,coupling_constraints,diet_fecal_compartments,diet,compress,write_lp,sample_label):
    model_out = model_dir+'%s.%s'%(sample_label,cobra_type.split('.')[1]) if not compress else model_dir+sample_label+'.xml.gz'
    lp_out = problem_dir+'%s.%s'%(sample_label,lp_type.split('.')[1])

    coverage_df = coverage_df.loc[coverage_df.sample_id==sample_label]
    original_id = coverage_df.original_id.unique()[0]
    pymgpipe_model = None

    if os.path.exists(model_out):
        pymgpipe_model = load_cobra_model(model_out)
        if write_lp:
            write_lp_problem(pymgpipe_model,out_file=lp_out,compress=compress,force=False)
    else:
        #pymgpipe_model = _build_com(sample_label=sample_label,tax=coverage_df,cutoff=1e-6,solver=solver)
        pymgpipe_model = build(
            name=sample_label,
            taxonomy=coverage_df,
            rel_threshold=1e-6,
            solver=solver,
            coupling_constraints=coupling_constraints,
            diet_fecal_compartments=diet_fecal_compartments
        )
        if diet is not None:
            personalized = load_dataframe(diet,return_empty=True) # try loading personlized diet
            if original_id in personalized.columns:
                diet = personalized[original_id].to_frame()
        
            add_diet_to_model(pymgpipe_model,diet)
        write_cobra_model(pymgpipe_model,model_out)
        if write_lp:
            write_lp_problem(pymgpipe_model,out_file=lp_out,compress=compress,force=True)
            
    del pymgpipe_model
    gc.collect()
    return model_out

def _format_coverage_file(coverage_file,taxa_dir,out_dir):
    existing_taxa_files = {t.split('/')[-1].split('.')[0]:taxa_dir+t for t in os.listdir(taxa_dir)}

    coverage = pd.read_csv(coverage_file,index_col=0,header=0)

    conversion_file_path = out_dir+'sample_label_conversion.csv'
    if not os.path.exists(conversion_file_path):
        sample_conversion_dict = {v:'mc'+str(i+1) for i,v in enumerate(coverage.columns)}
    else:
        sample_conversion_dict = pd.read_csv(conversion_file_path,index_col=0).iloc[:, 0].to_dict()
        if set(sample_conversion_dict.keys()) != set(coverage.columns):
            raise Exception('Provided label conversion file %s does not provide labels for all samples!'%conversion_file_path)
    
    conversion_t = {v:k for k,v in sample_conversion_dict.items()}
    coverage.rename(columns=sample_conversion_dict,inplace=True)
    
    sample_conversion_dict = pd.DataFrame({'conversion':sample_conversion_dict})
    sample_conversion_dict.to_csv(conversion_file_path)

    melted = pd.melt(coverage, value_vars=coverage.columns,ignore_index=False,var_name='sample_id',value_name='abundance',)
    melted['strain']=melted.index
    melted['file']=melted.strain.apply(lambda x: existing_taxa_files[x] if x in existing_taxa_files else None)
    melted['original_id'] = melted.sample_id.apply(lambda x: conversion_t[x])

    missing_taxa = melted[melted.file.isna()].index.unique()
    if len(missing_taxa)>0:
        melted.drop(missing_taxa,axis='index',inplace=True)
        print('Removed %s missing taxa from coverage file-' %len(missing_taxa))
        print(missing_taxa)
    
    melted = melted.loc[melted.abundance!=0]
    melted.abundance = melted.abundance.div(melted.groupby(['sample_id'])['abundance'].transform('sum'))

    melted.reset_index(inplace=True)
    melted.rename({'index':'id'},axis='columns',inplace=True)

    return melted

# def _build_com(sample_label, tax, cutoff, solver):
#     from micom.community import Community
    
#     com = Community(taxonomy=tax, progress=False, rel_threshold=cutoff, solver=solver,name=sample_label)
#     _add_pymgpipe_constraints(com=com,solver=solver)
#     com.solver.update()
#     return com

# def _add_pymgpipe_constraints(file=None,com=None,solver='gurobi'):
#     cobra_config = cobra.Configuration()
#     cobra_config.solver = solver

#     if com is None:
#         if file is None:
#             raise Exception('Need to pass in either file or model!')
#         try:
#             com = pickle.load(file)
#         except Exception as e:
#             return None
    
#     com.solver = solver
    
#     # RESETTING COEFS OF METAB EXCHANGE REACTIONS TO 1
#     exchange_reactions = [r for r in com.reactions if 'EX_' in r.id and 'biomass' not in r.id]
#     for react in exchange_reactions:
#         react.bounds = (-1000,1000)

#         metabs_to_remove = {metab[0]: metab[1] for metab in react.metabolites.items() if metab[1] != -1}
#         metabs_to_add = {metab[0]: 1.0 for metab in react.metabolites.items() if metab[1] != -1}
        
#         react.subtract_metabolites(metabs_to_remove)
#         react.add_metabolites(metabs_to_add)
    
#     abundances = com.abundances.to_dict()
#     biomass_reactions = _get_biomass_reactions(com)

#     biomass_and_biomass_exchange_reactions = [r for r in com.reactions if 'biomass' in r.id]
#     for reaction in biomass_and_biomass_exchange_reactions:
#         microbe = reaction.id.split('__')[1]
#         abundance = abundances[microbe]
#         reaction.bounds=(abundance,abundance)

#     com.solver.remove(com.solver.constraints.community_objective_equality)
#     com.solver.remove(com.solver.variables.community_objective)
    
#     expr = sum([r.flux_expression for r in biomass_reactions])
#     objective_constraint = com.problem.Constraint(name='objective_constraint',expression=expr,lb=0,ub=1000)
#     com.solver.add(objective_constraint)
#     com.objective = com.problem.Objective(expression=expr,direction='max')

#     _add_real_coupling_constraints(com)

#     return com

# def _add_coupling_constraints(com,u_const=0.01,C_const=400):
#     abundances = _get_model_abundances(com)
#     target_reactions = [r for r in com.reactions if not (('EX_' in r.id and '(e)' not in r.id and '[e]' not in r.id) or 'biomass' in r.id or 'Community' in r.id)]
#     for r in target_reactions:
#         split = ('___' if '___' in r.id else '__')
#         taxon = r.id.split(split)[1]
#         abundance = abundances[taxon]
#         lb = (-abundance*C_const) - u_const if r.lower_bound != 0 else 0
#         ub = (abundance*C_const) + u_const if r.upper_bound != 0 else 0

#         r.bounds = (lb,ub)

# def _add_real_coupling_constraints(com,u_const=0.01,C_const=400):
#     biomass_rxns= {r.community_id:r for r in com.reactions.query(re.compile('^biomass.*'))}

#     consts = []
#     target_reactions = [r for r in com.reactions if not (('EX_' in r.id and '(e)' not in r.id and '[e]' not in r.id) or 'biomass' in r.id or 'Community' in r.id)]
#     for r in target_reactions:
#         taxon = r.community_id
#         abundance = com.variables[biomass_rxns[taxon].id]
        
#         forward = r.forward_variable
#         reverse = r.reverse_variable

#         consts.append(com.solver.interface.Constraint(forward-(abundance*C_const),ub=u_const,name='%s_cp'%forward.name))
#         consts.append(com.solver.interface.Constraint(reverse-(abundance*C_const),ub=u_const,name='%s_cp'%reverse.name))
    
#     com.solver.add(consts)
#     com.solver.update()

# def _get_model_abundances(com):
#     abundances = {}
#     bms_rxns = _get_biomass_reactions(com)

#     for r in bms_rxns:
#         taxon = r.id.split('__')[1]
#         abundances[taxon] = r.lower_bound
#     return abundances

# def _get_biomass_reactions(com):
#     return [r for r in com.reactions if 'biomass' in r.id and 'EX' not in r.id]

def _is_valid_lp(file):
    with open(file, "rb") as fh:
        fh.seek(-1024, 2)
        last = fh.readlines()[-1].decode()
        if file.endswith('.lp'):
            return last.strip()=='End'
        elif file.endswith('.mps'):
            return last.strip()=='ENDATA'
        else:
            raise Exception('Unrecognized LP file at %s. Must be either .lp or .mps!'%file)

def _is_valid_sbml(file):
    with open(file, "rb") as fh:
        fh.seek(-1024, 2)
        last = fh.readlines()[-1].decode()
        return last.strip()=='</sbml>'

def _mute():
    sys.stdout = open(os.devnull, 'w')  

