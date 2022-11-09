import pandas as pd
import os

from .fva import regularFVA
from .optlang_util import *

def compute_nmpcs(
    models=[],
    out_file = None,
    reactions=None,
    regex=None,
    ex_only=True,
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    parallel=True,
    diet_fecal_compartments=False
):
    combined = pd.DataFrame()
    for s in models:
        m = load_model(path=s,solver=solver) if isinstance(s,str) else s
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
            metabs = res.index.str.split('[').str[0].drop_duplicates()
            df = {}
            for metab in metabs:
                fe = res.loc[metab+'[fe]']['max']
                d = res.loc[metab+'[d]']['min']
                
                df[metab]=d+fe
            nmpc = pd.DataFrame({m.name:df})
        else:
            nmpc = res['min']+res['max']
            nmpc.name = m.name

        combined = pd.concat([combined,nmpc],axis=1)

    if out_file is not None:
        combined.to_csv(out_file)        
    return combined