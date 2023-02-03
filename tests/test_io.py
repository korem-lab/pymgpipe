import tempfile
import optlang
import pytest
from pkg_resources import resource_filename
from pytest_check import check
from pymgpipe import *


def test_load_xml_model():
    m = load_model(pytest.sbml_model_path)
    assert isinstance(m, optlang.Model) and len(m.variables) == 12182


def test_load_mps_model():
    m = load_model(pytest.mps_problem_path)
    assert isinstance(m, optlang.Model) and len(m.variables) == 12182


def test_load_cobra_model():
    m = load_cobra_model(pytest.sbml_model_path)
    assert isinstance(m, cobra.Model) and len(m.reactions) == 6091


def test_write_mps_model(small_cobra_model):
    with tempfile.NamedTemporaryFile(delete=True, suffix=".mps") as tmp:
        write_lp_problem(small_cobra_model, out_file=tmp.name, compress=False)
        with check:
            assert os.path.exists(tmp.name)

        loaded = load_model(tmp.name)
        assert len(loaded.variables) == len(small_cobra_model.variables)


def test_write_lp_model(small_cobra_model):
    with tempfile.NamedTemporaryFile(delete=True, suffix=".lp") as tmp:
        write_lp_problem(small_cobra_model, out_file=tmp.name, compress=False)
        with check:
            assert os.path.exists(tmp.name)

        loaded = load_model(tmp.name)
        assert len(loaded.variables) == len(small_cobra_model.variables)


def test_write_cobra_model(small_cobra_model):
    with tempfile.NamedTemporaryFile(delete=True, suffix=".xml.gz") as tmp:
        write_cobra_model(small_cobra_model, file=tmp.name)
        with check:
            assert os.path.exists(tmp.name)

        loaded = load_cobra_model(tmp.name)
        assert isinstance(loaded, cobra.Model)
