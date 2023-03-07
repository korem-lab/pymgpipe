from pymgpipe import *

def test_remove_reverse_reactions(mini_optlang_model):
    num_reactions = len(mini_optlang_model.variables)
    some_var = mini_optlang_model.variables[100]
    reverse_var_id = get_reverse_id(some_var.name)
    assert reverse_var_id in mini_optlang_model.variables

    remove_reverse_vars(mini_optlang_model)
    new_reactions = len(mini_optlang_model.variables)

    assert new_reactions == num_reactions/2 and reverse_var_id not in mini_optlang_model.variables
