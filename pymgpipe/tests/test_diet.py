from pkg_resources import resource_filename
from pymgpipe import *


def test_add_diet_cobra(mini_cobra_model):
    added = add_diet_to_model(mini_cobra_model, "AverageEuropeanDiet")
    remove_diet(mini_cobra_model)
    assert len(added) == 10


def test_add_diet_optlang(mini_optlang_model):
    added = add_diet_to_model(mini_optlang_model, "AverageEuropeanDiet")
    remove_diet(mini_optlang_model)
    assert len(added) == 10


def test_remove_diet(mini_cobra_model):
    add_diet_to_model(mini_cobra_model, "AverageEuropeanDiet")
    remove_diet(mini_cobra_model)

    curr_diet = get_diet(mini_cobra_model)
    assert curr_diet.lb.unique()[0] == -1000 and curr_diet.ub.unique()[0] == 1000


def test_personlized_diet_force_uptake(mini_cobra_model):
    my_diet = [
        ["EX_pi(e)", 30.96],
        ["EX_h2o(e)", 35.59],
        ["EX_gln__L(e)", 44.78],
        ["EX_pyr(e)", 26.32],
    ]
    diet = pd.DataFrame(data=my_diet, columns=["Reaction", "flux"])
    diet.set_index(["Reaction"], inplace=True)

    remove_diet(mini_cobra_model)
    add_diet_to_model(mini_cobra_model, diet)

    expected_bounds = [
        ["Diet_EX_h2o[d]", -35.59, -28.472],
        ["Diet_EX_pyr[d]", -26.32, -21.056],
        ["Diet_EX_gln__L[d]", -44.78, -35.824],
        ["Diet_EX_pi[d]", -30.96, -24.768],
    ]
    expected = pd.DataFrame(data=expected_bounds, columns=["id", "lb", "ub"])
    expected.set_index("id", inplace=True)

    added_diet = get_diet(mini_cobra_model).set_index("id").loc[expected.index]

    remove_diet(mini_cobra_model)
    assert expected.round(3).equals(added_diet.loc[expected.index].round(3))


def test_personalized_diet_no_force_uptake(mini_cobra_model):
    my_diet = [
        ["EX_pi(e)", 30.96],
        ["EX_h2o(e)", 35.59],
        ["EX_gln__L(e)", 44.78],
        ["EX_pyr(e)", 26.32],
    ]
    diet = pd.DataFrame(data=my_diet, columns=["Reaction", "flux"])
    diet.set_index(["Reaction"], inplace=True)

    remove_diet(mini_cobra_model)
    add_diet_to_model(mini_cobra_model, diet, force_uptake=False)

    expected_bounds = [
        ["Diet_EX_h2o[d]", -35.59, 0],
        ["Diet_EX_pyr[d]", -26.32, 0],
        ["Diet_EX_gln__L[d]", -44.78, 0],
        ["Diet_EX_pi[d]", -30.96, 0],
    ]
    expected = pd.DataFrame(data=expected_bounds, columns=["id", "lb", "ub"])
    expected.set_index("id", inplace=True)

    added_diet = get_diet(mini_cobra_model).set_index("id").loc[expected.index]

    remove_diet(mini_cobra_model)
    assert expected.round(3).equals(added_diet.loc[expected.index].round(3))
