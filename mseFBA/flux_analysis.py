from pymgpipe.utils import get_reverse_id, get_reactions
from optlang.symbolics import Zero
import numpy as np
import logging

def get_mse_expression(model, flux_map):
    obj_expression = Zero

    flux_map = {k:v for k,v in flux_map.items() if k in model.variables and not np.isinf(v)}
    for f_id, flux in flux_map.items():
        forward_var = model.variables[f_id]
        reverse_var = model.variables[get_reverse_id(f_id)]
        net = forward_var-reverse_var

        squared_diff=(net-flux)**2
        obj_expression = obj_expression + squared_diff
    if obj_expression is Zero:
        logging.warning('No metabolites in mseFBA correlation objective, returning None')
    return obj_expression

def get_variance_expression(model, ids, weighted = True):
    vrs = [(v,model.variables[get_reverse_id(v.name)]) for v in get_reactions(model,reactions=ids)]
    if len(vrs) <= 1:
        logging.warning('No metabolites in variance objective, returning Zero')
        return Zero

    try:
        model.variables[0].primal
    except:
        model.solve()
    weight = 1/abs(np.mean([(f.primal-r.primal) for f,r in vrs]))
    weight = 1 if not weighted or np.isnan(weight) or np.isinf(weight) else weight

    net = [f-r for f,r in vrs]
    mean = np.sum(net)/len(net)
    obj_expr =  weight * np.sum([(v-mean)*(v-mean) for v in net])
    return obj_expr

def get_pfba_expression(model,ids):
    import itertools
    vrs = list(itertools.chain([(v,model.variables[get_reverse_id(v)]) for v in get_reactions(model,reactions=ids)]))
    if len(vrs)==0:
        logging.warning('No metabolites in pFBA objective, returning Zero')
        return Zero
    return np.sum(vrs)
