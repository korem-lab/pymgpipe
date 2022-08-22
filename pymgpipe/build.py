import sys
sys.path.insert(1, '../')

import pandas as pd
import os
from micom.community import Community
from micom import load_pickle
import cobra
from cobra.io import write_sbml_model
from pathlib import Path
import ctypes

from multiprocessing import Pool
from functools import partial
import gc
import tqdm
from cobra.io import read_sbml_model
from contextlib import redirect_stdout, contextmanager


def build_models(
    coverage_file,
    taxa_dir,
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    samples=None,
    parallelize=True
):
    gc.enable()

    if taxa_dir[-1] != '/':
        taxa_dir=taxa_dir+'/'

    Path('models').mkdir(exist_ok=True)
    Path('problems').mkdir(exist_ok=True)

    formatted_coverage_file=coverage_file.split('.csv')[0]+'_formatted.csv'
    if os.path.exists(formatted_coverage_file):
        formatted = pd.read_csv(formatted_coverage_file,index_col=0)
    else:
        formatted = _format_coverage_file(coverage_file,taxa_dir)
        formatted.to_csv(formatted_coverage_file)
    samples_to_run = formatted.sample_id.unique()
    taxa = formatted.strain.unique()

    print('Found coverage file with %s samples and %s unique taxa'%(len(samples_to_run),len(taxa)))
   
    if samples is not None:
        samples_to_run = samples if isinstance(samples,list) else [samples]
    
    threads = os.cpu_count() if threads == -1 else threads
    _func = partial(_build_single_model, formatted, solver)

    if parallelize:
        print('Building %s samples in parallel using %s threads...'%(len(samples_to_run),threads))
        p = Pool(threads, initializer=_mute)
        p.daemon = False

        built = list(tqdm.tqdm(p.imap(_func, samples_to_run),total=len(samples_to_run)))

        p.close()
        p.join()
    else:
        print('Building %s samples in series...'%(len(samples_to_run),threads))
        built = tqdm.tqdm(list(map(_func,samples_to_run)),total=len(samples_to_run))
    
    print('Finished building %s models and associated LP problems!'%len(built))    

def _build_single_model(coverage_df,solver,sample_label):
    model_out = 'models/%s.xml'%sample_label
    problem_out = 'problems/%s.mps'%sample_label
    coverage_df = coverage_df.loc[coverage_df.sample_id==sample_label]
    pymgpipe_model = None

    if os.path.exists(model_out) and _is_valid_sbml(model_out) and os.path.exists(problem_out) and _is_valid_lp(problem_out):
        return

    if os.path.exists(model_out) and _is_valid_sbml(model_out):
        pymgpipe_model = read_sbml_model(model_out)
        pymgpipe_model.solver.problem.write(problem_out)
        del pymgpipe_model
        gc.collect()
        return
    
    if not os.path.exists(model_out):
        pymgpipe_model = _build_com(sample_label=sample_label,tax=coverage_df,cutoff=1e-6,solver=solver)
        write_sbml_model(pymgpipe_model,model_out)
        pymgpipe_model.solver.problem.write(problem_out)

    del pymgpipe_model    
    gc.collect()
    return model_out

def _build_com(sample_label, tax, cutoff, solver):
    com = Community(taxonomy=tax, progress=False, rel_threshold=cutoff, solver=solver,name=sample_label)
    modified_com = _add_pymgpipe_constraints(com=com,solver=solver)
    return modified_com

def _add_pymgpipe_constraints(file=None,com=None,solver='gurobi'):
    cobra_config = cobra.Configuration()
    cobra_config.solver = solver

    if com is None:
        if file is None:
            raise Exception('Need to pass in either file or model!')
        try:
            com = load_pickle(file)
        except Exception as e:
            return None
    
    com.solver = solver
    
    # RESETTING COEFS OF METAB EXCHANGE REACTIONS TO 1
    exchange_reactions = [r for r in com.reactions if 'EX_' in r.id and 'biomass' not in r.id]
    for react in exchange_reactions:
        react.bounds = (-1000,1000)

        metabs_to_remove = {metab[0]: metab[1] for metab in react.metabolites.items() if metab[1] != -1}
        metabs_to_add = {metab[0]: 1.0 for metab in react.metabolites.items() if metab[1] != -1}
        
        react.subtract_metabolites(metabs_to_remove)
        react.add_metabolites(metabs_to_add)
    
    abundances = com.abundances.to_dict()
    biomass_reactions = _get_biomass_reactions(com)

    biomass_and_biomass_exchange_reactions = [r for r in com.reactions if 'biomass' in r.id]
    for reaction in biomass_and_biomass_exchange_reactions:
        microbe = reaction.id.split('__')[1]
        abundance = abundances[microbe]
        reaction.bounds=(abundance,abundance)

    com.solver.remove(com.solver.constraints.community_objective_equality)
    com.solver.remove(com.solver.variables.community_objective)
    
    expr = sum([r.flux_expression for r in biomass_reactions])
    objective_constraint = com.problem.Constraint(name='objective_constraint',expression=expr,lb=0,ub=1000)
    com.solver.add(objective_constraint)
    com.objective = com.problem.Objective(expression=expr,direction='max')

    _add_coupling_constraints(com)

    return com

def _format_coverage_file(coverage_file,taxa_dir):
    model_type = '.'+os.listdir(taxa_dir)[0].split('.')[1]
    coverage = pd.read_csv(coverage_file,index_col=0,header=0)
    sample_conversion_dict = {v:'mc'+str(i+1) for i,v in enumerate(coverage.columns)}
    coverage.rename(columns=sample_conversion_dict,inplace=True)
    
    sample_conversion_dict = pd.DataFrame({'conversion':sample_conversion_dict})
    sample_conversion_dict.to_csv('sample_label_conversion.csv')

    melted = pd.melt(coverage, value_vars=coverage.columns,ignore_index=False,var_name='sample_id',value_name='abundance',)
    melted['strain']=melted.index
    melted['file']=taxa_dir+melted.strain+model_type

    missing_taxa = set([f.split('.')[0] for f in melted.file if not os.path.exists(f)])
    if len(missing_taxa)>0:
        melted.drop(missing_taxa,axis='index',inplace=True)
        print('!!! Removed %s missing taxa from coverage file\n' %len(missing_taxa))
        print(missing_taxa)
    
    melted = melted.loc[melted.abundance!=0]

    melted.reset_index(inplace=True)
    melted.rename({'index':'id'},axis='columns',inplace=True)

    return melted

def _add_coupling_constraints(com,u_const=0.01,C_const=400):
    abundances = _get_model_abundances(com)
    target_reactions = [r for r in com.reactions if not (('EX_' in r.id and '(e)' not in r.id and '[e]' not in r.id) or 'biomass' in r.id or 'Community' in r.id)]
    for r in target_reactions:
        split = ('___' if '___' in r.id else '__')
        taxon = r.id.split(split)[1]
        abundance = abundances[taxon]
        lb = (-abundance*C_const) - u_const if r.lower_bound != 0 else 0
        ub = (abundance*C_const) + u_const if r.upper_bound != 0 else 0

        r.bounds = (lb,ub)

def _get_model_abundances(com):
    abundances = {}
    bms_rxns = _get_biomass_reactions(com)

    for r in bms_rxns:
        taxon = r.id.split('__')[1]
        abundances[taxon] = r.lower_bound
    return abundances

def _get_biomass_reactions(com):
    return [r for r in com.reactions if 'biomass' in r.id and 'EX' not in r.id]

def _is_valid_lp(file):
    with open(file, "rb") as fh:
        fh.seek(-1024, 2)
        last = fh.readlines()[-1].decode()
        return last.strip()=='ENDATA'

def _is_valid_sbml(file):
    with open(file, "rb") as fh:
        fh.seek(-1024, 2)
        last = fh.readlines()[-1].decode()
        return last.strip()=='</sbml>'

def _mute():
    sys.stdout = open(os.devnull, 'w')  