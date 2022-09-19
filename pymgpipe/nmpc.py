import pandas as pd
import os

from .fva import regularFVA
from .optlang_util import load_model

def compute_nmpcs(
    models=[],
    out_file = None,
    reactions=None,
    regex=None,
    ex_only=True,
    solver='gurobi',
    threads=int(os.cpu_count()/2),
    parallel=True
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
        nmpc = res['min']+res['max']
        nmpc.name = m.name

        combined = pd.concat([combined,nmpc],axis=1)
    if out_file is not None:
        combined.to_csv(out_file)
    return combined
