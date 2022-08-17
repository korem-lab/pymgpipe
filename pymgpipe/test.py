from .optlang_util import load_model, solve_model
from .fva import regularFVA

def test_pymgpipe(solver='gurobi'):
    loaded_model = load_model('sample_model.mps',solver=solver)
    print('Loaded model...')
    solution = solve_model(model=loaded_model)
    print('Solved model...')
    fva_sol = regularFVA(path='sample_model.mps',ex_only=True,solver=solver)
    print('Performed FVA...\n')

    print(f'Finished testing! All good to go!')
