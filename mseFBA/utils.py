import os
import pandas as pd
import logging
from pymgpipe import load_dataframe

def evaluate_results(to_compare,metabolomics,axis=1):
    to_compare = load_dataframe(to_compare)
    metabolomics = load_dataframe(metabolomics)

    if len(metabolomics.columns) != len(to_compare.columns):
        logging.warning('Comparing two dataframes with unequal number of samples- %s vs %s'%(len(to_compare.columns),len(metabolomics.columns)))

    sp_r = to_compare.corrwith(metabolomics,method='spearman',axis=axis)
    pr_r = to_compare.corrwith(metabolomics,method='pearson',axis=axis)
    combined_dataframe = pd.concat([sp_r,pr_r],axis=1)
    combined_dataframe.columns=['spearman','pearson']

    print('Avg Spearman- '+str(sp_r.mean()))
    print('Avg Pearson- '+str(pr_r.mean()))

    return combined_dataframe

def format_nmpc_from_files(fva_dir='fva/',out_file='nmpc_sol.csv',write_to_file=True):
    fva_files = os.listdir(fva_dir)
    print('Computing NMPCs using fva results for %s samples...'%len(fva_files))
    
    nmpcs = {}
    for f in fva_files:
        if 'min_max' in f:
            continue
        s = load_dataframe(fva_dir+f)
        label = f.split('.csv')[0]

        nmpcs[label]=(s['min']+s['max']).to_dict()

    nmpcs = pd.DataFrame(nmpcs)
    nmpcs.sort_index(axis=1,inplace=True)
    nmpcs.sort_index(axis=0,inplace=True)
    if write_to_file:
        nmpcs.to_csv(out_file)
    return nmpcs

def get_objective_value(m):
    import optlang
    if m.interface is optlang.cplex_interface:
        return m.objective.expression.as_coefficients_dict()[1] + m.problem.solution.get_objective_value()
    elif m.interface is optlang.gurobi_interface:
        return m.objective.expression.as_coefficients_dict()[1] + m.problem.getAttr("ObjVal")
    else:
        raise Exception('Unrecognized solver- %s'%m.interface)