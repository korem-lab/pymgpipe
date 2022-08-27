import os
import pandas as pd

def evaluate_results(to_compare,metabolomics):
    to_compare = load_dataframe(to_compare)
    metabolomics = load_dataframe(metabolomics)

    sp_r = to_compare.corrwith(metabolomics,method='spearman',axis=1)
    pr_r = to_compare.corrwith(metabolomics,method='pearson',axis=1)
    combined_dataframe = pd.concat([sp_r,pr_r],axis=1)
    combined_dataframe.columns=['spearman','pearson']

    print('Avg Spearman- '+str(sp_r.mean()))
    print('Avg Pearson- '+str(pr_r.mean()))

    return combined_dataframe

def compute_nmpcs(fva_dir='fva/',out_file='nmpc_sol.csv',write_to_file=True):
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

def load_dataframe(m, return_empty=False):
    if isinstance(m,str):
        if not os.path.exists(m):
            if return_empty:
                return pd.DataFrame() 
            else:
                raise Exception('Tried to load dataframe from path that does not exist- %s'%m)
        
        return pd.read_csv(m,index_col=0)
    elif isinstance(m, pd.DataFrame):
        return m
    else:
        raise Exception('_load_dataframe can only take a string or dataframe, received %s'%type(m))
