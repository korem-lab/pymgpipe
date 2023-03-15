import os
import pandas as pd
import csv
import pathlib
import contextlib
from .io import load_model, suppress_stdout
from .utils import get_reactions, Constants

class VFFVA(object):
    def __init__(self):
        self._path = None

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        self._path = value

vffva_config = VFFVA()

def veryFastFVA(
        nCores, 
        nThreads,
        path, 
        model=None, 
        ex_only=True,
        scaling=0, 
        memAff='none', 
        schedule='dynamic', 
        nChunk=50, 
        optPerc=90, 
        reactions=None,
        regex=None,
        threshold=1e-5):
    '''
    VFFVA performs Very Fast Flux Variability Analysis (VFFVA). VFFVA is a parallel implementation of FVA that
    allows dynamically assigning reactions to each worker depending on their computational load
    Guebila, Marouen Ben. "Dynamic load balancing enables large-scale flux variability analysis." bioRxiv (2018): 440701.

    USAGE:
    minFlux,maxFlux=VFFVA(nCores, nThreads, model, scaling, memAff, schedule, nChunk, optPerc, ex)

    :param nCores:   Number of non-shared memory cores/machines.
    :param nThreads: Number of shared memory threads in each core/machine.
    :param model:    .mps format: path to metabolic model in .mps format.
    :param scaling:  CPLEX parameter. It corresponds to SCAIND parameter (Default = 0).
                     -1: no scaling; 0: equilibration scaling; 1: more aggressive scaling.
                     more information here: https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/ScaInd.html.
    :param memAff:   'none', 'core', or 'socket'. (Default = 'none'). This an OpenMPI parameter, more
                     information here: https://www.open-mpi.org/faq/?category=tuning#using-paffinity-v1.4.
    :param schedule: 'dynamic', 'static', or 'guided'. (Default = 'dynamic'). This is an OpenMP parameter, more
                     information here: https://software.intel.com/en-us/articles/openmp-loop-scheduling
    :param nChunk:   Number of reactions in each chunk (Default = 50). This is an OpenMP parameter, more
                     information here: https://software.intel.com/en-us/articles/openmp-loop-scheduling
    :param optPerc:  Percentage of the optimal objective used in FVA. Integer between 0 and 100. For example, when set to 90
                     FVA will be computed on 90% of the optimal objective.
    :param ex:      0-based indices of reactions to optimize as a numpy array.  (Default, all reactions)
    :return:
           minFlux:          (n,1) vector of minimal flux values for each reaction.
           maxFlux:          (n,1) vector of maximal flux values for each reaction.
    '''
    status = os.system('mpirun --version')
    if status != 0:
        raise ValueError(['MPI and/or CPLEX nont installed, please follow the install guide'
               'or use the quick install script'])

    # Set schedule and chunk size parameters
    os.environ["OMP_SCHEDUELE"] = schedule+str(nChunk)

    model = load_model(path if model is None else model)
    var_dict = {i:v.name for i,v in enumerate(model.variables)}
    var_dict_inv = {v:k for k,v in var_dict.items()}

    ex_indices = [var_dict_inv[r.name] for r in get_reactions(model, reactions = reactions, regex = Constants.EX_REGEX if ex_only and regex is None else regex)]
    print('Running on %s reactions!'%len(ex_indices))
    
    # Set reactions to optimize
    rxns_file = 'rxns_%s.txt'%model.name
    with open(rxns_file, "w") as f:
        for num in ex_indices:
            f.write(str(num) + "\n")   
  
    assert vffva_config.path is not None, 'Please set value of `vffva_config.path` to location of VFFVA executable'
    status = os.system('mpirun -np ' + str(nCores) + ' --bind-to ' + str(memAff) + ' -x OMP_NUM_THREADS=' + str(nThreads) +
        f' {vffva_config.path}/lib/veryfastFVA ' + path + ' ' + str(optPerc) + ' ' + str(scaling) + ' ' + rxns_file)

    # Fetch results
    resultFile = path[:-4] + 'output.csv'
    if not os.path.exists(resultFile):
        raise Exception('Ran into issue when running VFFVA, could not find results file...')
    results = pd.read_csv(resultFile,header=0)

    os.system('rm '+resultFile)
    os.system('rm '+rxns_file)
    results.rename({i:var_dict[x] for i,x in enumerate(ex_indices) if i in results.index},axis=0,inplace=True)

    results.index.rename('id',inplace=True)
    if threshold is not None:
        results[abs(results)<threshold]=0
    
    results.sort_index(inplace=True)
    results.columns=['min','max']
    return results
