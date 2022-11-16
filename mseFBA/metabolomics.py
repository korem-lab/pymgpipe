from .utils import *
import numpy as np
import logging

def process_metabolomics(metabolomics,fva_dir='fva/',scale=True,map_labels=True,conversion_file='sample_label_conversion.csv',out_file=None):
    metabolomics_df = load_dataframe(metabolomics)

    metabolomics_df = _map_sample_labels(metabolomics_df,conversion_file) if map_labels else metabolomics_df
    metabolomics_df = _scale_metabolomics(metabolomics_df,fva_dir) if scale else metabolomics_df

    print('\n------------------------------------------------------------------')
    print('Using processed metabolomics with %s metabolites and %s samples!'%(len(metabolomics_df.index),len(metabolomics_df.columns)))
    print('------------------------------------------------------------------')
    if out_file is not None:
        metabolomics_df.to_csv(out_file)
    return metabolomics_df

def _map_sample_labels(metabolomics, conversion_file):
    metabolomics_df = load_dataframe(metabolomics)
    conversion = load_dataframe(conversion_file).conversion.to_dict()

    samples = list(set(metabolomics_df.columns).intersection(set(conversion.keys())))
    print('\nMapping sample labels for %s matched samples...'%len(samples))

    metabolomics_df = metabolomics_df[samples]
    return metabolomics_df.rename(conversion, axis='columns')

def _scale_metabolomics(metabolomics,fva_dir='fva/'):
    raw = load_dataframe(metabolomics)

    if not os.path.exists(fva_dir):
        raise Exception('Could not find fva directory at %s'%fva_dir)

    min_max = load_dataframe(fva_dir+'min_max.csv',return_empty=True)
    if min_max.empty:
        fva_dfs = [load_dataframe(fva_dir+f) for f in os.listdir(fva_dir)]
    
        if len(fva_dfs)==0:
            raise Exception('No FVA results found')

        print('\nScaling metabolomics using FVA results for %s out of %s total samples!'%(len(fva_dfs),len(raw.columns)))
        if len(fva_dfs)/len(raw.columns) < 0.5:
            logging.warning(f'Using FVA results for less than 50% of all samples, results might not be as good as they could be!\n')

        mins = pd.concat([df['min'] for df in fva_dfs],axis=1)
        maxs = pd.concat([df['max'] for df in fva_dfs],axis=1)

        min_max = pd.concat([mins.min(axis=1,numeric_only=True).rename('min'),maxs.max(axis=1,numeric_only=True).rename('max')],axis=1)
        min_max.to_csv(fva_dir+'min_max.csv')
    else:
        print('Using existing min-max FVA results! If you want to re-calculate ranges, delete %s'%(fva_dir+'min_max.csv'))
    
    min_max= min_max.T.to_dict()
    
    missing_metabs = [m for m in list(raw.index) if m not in min_max]
    if len(missing_metabs)>0:
        print('Skipping %s metabolites due to them being missing from models, these metabolites will not be included in formatted metabolomics...'%len(missing_metabs)) 
        raw.drop(missing_metabs,inplace=True)

    def _scale_row(x):
        metab = x.name
        log_row = np.log10(x)

        metab_min = min_max[metab]['min']
        metab_max = min_max[metab]['max']
        scaled = ((metab_max-metab_min)*(log_row-log_row.min())/(log_row.max()-log_row.min()))+metab_min
        return scaled
    
    scaled = raw.apply(_scale_row,axis=1)
    scaled.sort_index(inplace=True)

    if len(scaled.index) == 0:
        raise Exception('Scaled metabolomics is empty. Please check FVA files to make sure everything is correct!')

    return scaled

class transformations:
    def none(x):
        return x
    def log10(x):
        return np.log10(x)
    def rzsc(x, q = 0.05):
        return (x-x.median()) / (x[(x >= x.quantile(q)) & (x <= x.quantile(1-q))]).std()
    def abundance(x):
        log_row = np.log10(x)
        return log_row/sum(log_row.dropna())
    def minmax(x,min_val=-1000,max_val=1000):
        log_row = np.log10(x)
        scaled = ((max_val-min_val)*(log_row-log_row.min())/(log_row.max()-log_row.min()))+min_val
        return scaled

def transform_metabolomics(metabolomics,func=transformations.none,*args):
    raw = load_dataframe(metabolomics)   

    print('Transforming metabolomics with %s samples and %s metabolites...'%(len(raw.columns),len(raw.index)))
    if 'transformations.' in transformations.rzsc.__qualname__:
        print('Using pre-built transformation function- %s'%str(func.__name__))

    axis = 0 if func is transformations.abundance else 1 
    scaled = raw.apply(func,args=args,axis=axis)

    return scaled

