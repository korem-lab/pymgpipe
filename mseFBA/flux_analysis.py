from pymgpipe.optlang_util import get_reverse_id, get_reactions
from optlang.symbolics import Zero
import numpy as np
import logging

def get_mse_expression(model, flux_map):
    obj_expression = Zero

    flux_map = {k:v for k,v in flux_map.items() if k in model.variables}
    for f_id, flux in flux_map.items():
        forward_var = model.variables[f_id]
        reverse_var = model.variables[get_reverse_id(f_id)]
        net = forward_var-reverse_var

        squared_diff=(net-flux)**2
        obj_expression = obj_expression + squared_diff
    if obj_expression is None:
        logging.warning('No metabolites in mseFBA correlation objective, returning None')
    return obj_expression

def get_variance_expression(model, ids):
    vrs = [v-model.variables[get_reverse_id(v.name)] for v in get_reactions(model,reactions=ids)]
    if len(vrs) <= 1:
        logging.warning('No metabolites in variance objective, returning Zero')
        return Zero

    mean = np.sum(vrs)/len(vrs)
    obj_expr = np.sum([(v-mean)*(v-mean) for v in vrs])
    return obj_expr

def get_pfba_expression(model,ids):
    import itertools
    vrs = list(itertools.chain([(v,model.variables[get_reverse_id(v)]) for v in get_reactions(model,reactions=ids)]))
    if len(vrs)==0:
        logging.warning('No metabolites in pFBA objective, returning Zero')
        return Zero
    return np.sum(vrs)
