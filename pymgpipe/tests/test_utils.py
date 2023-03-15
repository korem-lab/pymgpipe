from pymgpipe import *

def test_remove_reverse_reactions(mini_optlang_model):
    num_reactions = len(mini_optlang_model.variables)
    some_var = mini_optlang_model.variables[100]
    reverse_var_id = get_reverse_id(some_var.name)
    assert reverse_var_id in mini_optlang_model.variables

    remove_reverse_vars(mini_optlang_model,hard_remove=False)
    assert len(mini_optlang_model.variables) == num_reactions and \
        mini_optlang_model.variables[reverse_var_id].lb == 0 and mini_optlang_model.variables[reverse_var_id].ub == 0

    remove_reverse_vars(mini_optlang_model,hard_remove=True)
    assert len(mini_optlang_model.variables) == num_reactions/2 and reverse_var_id not in mini_optlang_model.variables


