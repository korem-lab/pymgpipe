from pkg_resources import resource_filename
from pymgpipe import *


def test_add_diet(small_cobra_model):
    added = add_diet_to_model(small_cobra_model, "AverageEuropeanDiet")
    get_adapted_diet
    assert not added.empty


def test_remove_diet(small_cobra_model):
    add_diet_to_model(small_cobra_model, "AverageEuropeanDiet")
    remove_diet(small_cobra_model)

    curr_diet = get_diet(small_cobra_model)
    assert curr_diet.lb.unique()[0] == -1000 and curr_diet.ub.unique()[0] == 1000


def test_personlized_diet(small_cobra_model):
    my_diet = [
        ["EX_his_L(e)", 30.96],
        ["EX_etoh(e)", 35.59],
        ["EX_h2o(e)", 44.78],
        ["EX_fe2(e)", 26.32],
    ]
    diet = pd.DataFrame(data=my_diet, columns=["Reaction", "small_sample_model"])
    diet.set_index(["Reaction"], inplace=True)
    remove_diet(small_cobra_model)
    added = add_diet_to_model(small_cobra_model, diet)

    expected_bounds = [
        ["Diet_EX_etoh[d]", -35.59, -28.47],
        ["Diet_EX_fe2[d]", -26.32, -21.056],
        ["Diet_EX_h2o[d]", -44.78, -35.824],
        ["Diet_EX_his_L[d]", -30.96, -24.77],
    ]
    expected = pd.DataFrame(data=expected_bounds, columns=["id", "lb", "ub"])
    expected.set_index("id", inplace=True)

    assert added.loc[added.ub != 0].set_index("id").round(2).equals(expected.round(2))
