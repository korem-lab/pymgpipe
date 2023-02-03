import pytest
from pymgpipe import load_cobra_model, load_model
from pkg_resources import resource_filename


def pytest_configure():
    pytest.sbml_model_path = resource_filename(
        "pymgpipe", "resources/models/small_sample_model.xml.gz"
    )
    pytest.mps_problem_path = resource_filename(
        "pymgpipe", "resources/problems/small_sample_model.mps.gz"
    )


@pytest.fixture(scope="session")
def small_cobra_model():
    m = load_cobra_model(pytest.sbml_model_path)
    return m


@pytest.fixture(scope="session")
def small_optlang_model():
    m = load_model(pytest.mps_problem_path)
    return m
