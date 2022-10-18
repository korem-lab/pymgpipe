from json import load
from .optlang_util import load_model, solve_model, get_reactions, Constants
from .fva import regularFVA
import pkg_resources

def test_pymgpipe(solver='gurobi'):
    print('Loading model...')
    loaded_model = get_single_sample_model(solver=solver)

    print('Fetching exchange reactions...')
    ex_reactions = get_reactions(model=loaded_model,regex=Constants.EX_REGEX)
    if len(ex_reactions)==0:
        raise Exception('Error fetching exchange reactions from model!')

    print('Solving model...')
    solution = solve_model(model=loaded_model,ex_only=True)
    if len(solution.index) != len(ex_reactions):
        raise Exception('Reguar FBA solution did not have expected number of metabolites!')

    print('Performing FVA on a small subset of reactions...')
    fva_sol = regularFVA(model=loaded_model,reactions=ex_reactions[:20],write_to_file=False)
    print(fva_sol)
    if len(fva_sol.index) != 20:
        raise Exception('FVA solution did not have expected number of metabolites!')

    print(f'Finished testing! All good to go!')

def get_single_sample_model(solver='gurobi'):
    resource_path = pkg_resources.resource_filename('pymgpipe','samples/single_sample_model.mps')
    loaded_model = load_model(resource_path,solver=solver)
    return loaded_model

def get_multi_sample_model(solver='gurobi'):
    resource_path = pkg_resources.resource_filename('pymgpipe','samples/multi_sample_model.mps')
    loaded_model = load_model(resource_path,solver=solver)
    return loaded_model