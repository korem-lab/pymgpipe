from pkg_resources import resource_filename
from pymgpipe import *


def test_add_diet():
    m = load_cobra_model(
        resource_filename("pymgpipe", "resources/models/small_sample_model.xml.gz")
    )
    added = add_diet_to_model(m, "AverageEuropeanDiet")

    assert not added.empty


def test_remove_diet():
    m = load_cobra_model(
        resource_filename("pymgpipe", "resources/models/small_sample_model.xml.gz")
    )
    add_diet_to_model(m, "AverageEuropeanDiet")
    remove_diet(m)

    curr_diet = get_diet(m)

    assert curr_diet.lb.unique()[0] == -1000 and curr_diet.ub.unique()[0] == 1000
