from .optlang_util import load_model, solve_model, _get_exchange_reactions
from .fva import regularFVA
import pkg_resources

def test_pymgpipe(solver='gurobi'):
    resource_path = pkg_resources.resource_filename('pymgpipe','sample_model.mps')

    print('Loading model...')
    loaded_model = load_model(resource_path,solver=solver)

    print('Fetching exchange reactions...')
    ex_reactions = _get_exchange_reactions(model=loaded_model)
    if len(ex_reactions)==0:
        raise Exception('Error fetching exchange reactions from model!')

    print('Solving model...')
    solution = solve_model(model=loaded_model,ex_only=True)
    if len(solution.index) != len(ex_reactions):
        raise Exception('Reguar FBA solution did not have expected number of metabolites!')

    print('Performing FVA on a small subset of reactions...')
    fva_sol = regularFVA(path=resource_path,ex_only=True,solver=solver,specific_reactions=ex_reactions[:20],write_to_file=False)
    if len(fva_sol.index) != 20:
        raise Exception('FVA solution did not have expected number of metabolites!')

    print(f'Finished testing! All good to go!')
