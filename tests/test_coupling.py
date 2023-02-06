import re
from pymgpipe import add_coupling_constraints, remove_coupling_constraints
from pytest_check import check


def test_add_coupling_optlang(mini_optlang_model):
    add_coupling_constraints(mini_optlang_model)

    added = [k for k in mini_optlang_model.constraints if re.match(".*_cp$", k.name)]
    assert len(added) > 0


def test_add_coupling_cobra(mini_cobra_model):
    add_coupling_constraints(mini_cobra_model)

    added = [k for k in mini_cobra_model.constraints if re.match(".*_cp$", k.name)]
    assert len(added) > 0


def test_remove_coupling_optlang(mini_optlang_model):
    remove_coupling_constraints(mini_optlang_model)

    before = len(mini_optlang_model.constraints)
    add_coupling_constraints(mini_optlang_model)

    with check:
        assert len(mini_optlang_model.constraints) > before

    remove_coupling_constraints(mini_optlang_model)
    assert len(mini_optlang_model.constraints) == before


def test_remove_coupling_cobra(mini_cobra_model):
    remove_coupling_constraints(mini_cobra_model)

    before = len(mini_cobra_model.constraints)
    add_coupling_constraints(mini_cobra_model)

    with check:
        assert len(mini_cobra_model.constraints) > before

    remove_coupling_constraints(mini_cobra_model)
    assert len(mini_cobra_model.constraints) == before
