import tempfile
import optlang
import pytest
from pkg_resources import resource_filename
from pytest_check import check
from pymgpipe import *
import random


lp_problem_ext = [".lp", ".lp.gz", ".mps", ".mps.gz"]
cobra_problem_ext = [".json", ".mat", ".xml", ".xml.gz"]


def test_read_cobra():
    for ext in cobra_problem_ext:
        m = load_cobra_model(pytest.resource_models_dir + "mini_model" + ext)
        with check:
            assert len(m.variables) == 926 and len(m.constraints) == 354


def test_read_lp():
    for ext in lp_problem_ext:
        m = load_model(pytest.resource_problems_dir + "mini_model" + ext)
        with check:
            assert len(m.variables) == 926 and len(m.constraints) == 354


def test_write_lp(mini_cobra_model):
    for s in lp_problem_ext:
        with tempfile.NamedTemporaryFile(delete=True, suffix=s) as tmp:
            write_lp_problem(mini_cobra_model, out_file=tmp.name, compress=False)
            with check:
                assert os.path.exists(tmp.name)

            loaded = load_model(tmp.name)
            assert len(loaded.variables) == len(mini_cobra_model.variables)


def test_write_cobra(mini_cobra_model):
    for s in cobra_problem_ext:
        with tempfile.NamedTemporaryFile(delete=True, suffix=s) as tmp:
            ex1 = mini_cobra_model.reactions.get_by_id('EX_h2o[fe]')
            ex2 = mini_cobra_model.reactions.get_by_id('EX_co2[fe]')
            ex2.bounds = (-5000,5000)

            write_cobra_model(mini_cobra_model, tmp.name)

            with check:
                assert os.path.exists(tmp.name)

            loaded = load_cobra_model(tmp.name)
            ex1_new = loaded.reactions.get_by_id('EX_h2o[fe]')
            ex2_new = loaded.reactions.get_by_id('EX_co2[fe]')

            assert ex1.bounds == ex1_new.bounds and ex2.bounds == ex2_new.bounds
            assert len(loaded.variables) == len(mini_cobra_model.variables)
