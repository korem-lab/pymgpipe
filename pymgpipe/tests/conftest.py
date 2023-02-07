import pytest
from pymgpipe import load_cobra_model, load_model
from pkg_resources import resource_filename


def pytest_configure():
    pytest.resource_models_dir = resource_filename("pymgpipe", "resources/models/")
    pytest.resource_problems_dir = resource_filename("pymgpipe", "resources/problems/")


@pytest.fixture(scope="session")
def mini_cobra_model():
    m = load_cobra_model(pytest.resource_models_dir + "mini_model.xml")
    return m


@pytest.fixture(scope="session")
def mini_optlang_model():
    m = load_model(pytest.resource_problems_dir + "mini_model.mps")
    return m
