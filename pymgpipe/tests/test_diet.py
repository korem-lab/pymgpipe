from pkg_resources import resource_filename
from pymgpipe import *


def test_add_diet(mini_cobra_model):
    added = add_diet_to_model(mini_cobra_model, "AverageEuropeanDiet")
    get_adapted_diet
    assert not added.empty


def test_remove_diet(mini_cobra_model):
    add_diet_to_model(mini_cobra_model, "AverageEuropeanDiet")
    remove_diet(mini_cobra_model)

    curr_diet = get_diet(mini_cobra_model)
    assert curr_diet.lb.unique()[0] == -1000 and curr_diet.ub.unique()[0] == 1000


def test_personlized_diet(mini_cobra_model):
    my_diet = [
        ["EX_pi(e)", 30.96],
        ["EX_h2o(e)", 35.59],
        ["EX_gln__L(e)", 44.78],
        ["EX_pyr(e)", 26.32],
    ]
    diet = pd.DataFrame(data=my_diet, columns=["Reaction", "flux"])
    diet.set_index(["Reaction"], inplace=True)
    remove_diet(mini_cobra_model)
    added = add_diet_to_model(mini_cobra_model, diet)

    expected_bounds = [
        ["Diet_EX_h2o[d]", -35.59, -28.472],
        ["Diet_EX_pyr[d]", -26.32, -21.056],
        ["Diet_EX_gln__L[d]", -44.78, -35.824],
        ["Diet_EX_pi[d]", -30.96, -24.768],
    ]
    expected = pd.DataFrame(data=expected_bounds, columns=["id", "lb", "ub"])
    expected.set_index("id", inplace=True)

    assert (
        expected.round(2)
        .sort_index()
        .equals(added.loc[added.ub != 0].set_index("id").round(2).sort_index())
    )
