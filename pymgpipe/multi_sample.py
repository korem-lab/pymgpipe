from .optlang_util import *
import tqdm
import sys
from multiprocessing import Pool
from functools import partial
from .optlang_util import _get_reverse_id
import optlang
import os
import numpy as np

def l_model(m_path,solver='gurobi'):
    m = load_model(m_path,solver=solver)
    return (m,'_'+m_path.split('/')[-1].split('.')[0])

def combine_samples(
    models,
    out_file=None,
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    parallel=False
    ):
    
    if isinstance(models,str):
        models = [models+m for m in os.listdir(models)]
    combined_model = optlang.Model()
    obj_expr = None

    threads = os.cpu_count()-1 if threads == -1 else threads
    threads = min(threads,len(models))

    print('Combining %s models using %s threads...'%(len(models),threads))
    sys.stdout = open(os.devnull, 'w')  

    _func = partial(l_model,solver=solver)
    if parallel:
        p = Pool(processes=threads)
        res = p.imap(_func, models)
    else:
        res = map(_func, models)

    for res in tqdm.tqdm(res,total=len(models)):
        curr,label=res
        label = label.split('#')[0] #sometimes a hashtag is added? weird
        for v in curr.variables:
            v.name = v.name + label
        for c in curr.constraints:
            c.name = c.name + label
        obj_expr = curr.objective.expression if obj_expr is None else obj_expr + curr.objective.expression

        combined_model._add_variables(curr.variables)
        combined_model._add_constraints(curr.constraints)
    combined_model.objective = optlang.Objective(obj_expr)
    # combined_model.update()
    if out_file is not None:
        combined_model.problem.write(out_file)
    sys.stdout = sys.__stdout__

    return combined_model


def solve_multi_sample_model(
    model,
    regex=None,
    reactions=None,
    solver='gurobi',
    verbosity=0,
    presolve=True,
    method='auto',
    flux_threshold=1e-5,
    ex_only=True):
    if isinstance(model,str):
        model = load_model(model, solver)
    model.configuration.verbosity=verbosity
    model.configuration.presolve=presolve
    model.configuration.lp_method=method
    model.configuration.qp_method=method

    model.update()
    
    try:
        model.optimize()
    except Exception as e:
        raise Exception(f"Error when optimizing model, make sure model is valid- ${e}")

    if model.status == 'infeasible':
        raise InfeasibleModelException('%s is infeasible!'%model.name)

    regex = '^EX_.*_m_.*$' if ex_only else regex
    fluxes = _get_fluxes_from_model(model,threshold=flux_threshold,regex=regex,reactions=reactions)
    return pd.DataFrame(fluxes)

def _get_fluxes_from_model(model,reactions=None,regex=None,threshold=1e-5):
    fluxes = {}
    rxns = get_reactions(model,reactions,regex)
    for forward in rxns:
        sample_id = 'mc'+forward.name.split('_mc')[1]
        r_id = _get_reverse_id(forward.name, multi_sample=True)
        if r_id not in model.variables:
            continue
        reverse = model.variables[r_id]

        flux = float(forward.primal-reverse.primal)
        flux = 0 if flux == -0.0 else flux
        flux = flux if abs(flux)>threshold else 0
        if sample_id not in fluxes:
            fluxes[sample_id] = {}
        fluxes[sample_id][forward.name.split('_mc')[0]]=flux
    return fluxes


def compute_multi_sample_nmpcs(model,reactions=None,ex_only=True):
    comb = pd.DataFrame()
    if reactions is None and ex_only is True:
        reactions = list(set([m.name.split('_mc')[0] for m in get_reactions(model,regex='^EX_.*_m_.*$')]))

    print('Performing FVA on %s reactions...'%len(reactions))
    for metab in tqdm.tqdm(reactions):
        metab_rxns = get_reactions(model,regex='%s_mc.*$'%metab)
        net_rxns = [f-model.variables[_get_reverse_id(f.name,multi_sample=True)] for f in metab_rxns]
        s_obj = np.sum(net_rxns)
        
        model.objective = optlang.Objective(s_obj,direction='min')
        min_sol = solve_multi_sample_model(model,reactions=[m.name for m in metab_rxns])
        
        model.objective = optlang.Objective(s_obj,direction='max')
        max_sol = solve_multi_sample_model(model,reactions=[m.name for m in metab_rxns])

        comb = pd.concat([comb,min_sol.add(max_sol)],axis=0)
    return comb


