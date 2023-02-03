from pkg_resources import resource_filename
import optlang
from pymgpipe import *


def test_load_xml_model():
    m = load_model(
        resource_filename("pymgpipe", "resources/models/small_sample_model.xml.gz")
    )
    assert isinstance(m, optlang.Model) and len(m.variables) == 12182


def test_load_mps_model():
    m = load_model(
        resource_filename("pymgpipe", "resources/problems/small_sample_model.mps.gz")
    )
    assert isinstance(m, optlang.Model) and len(m.variables) == 12182


def test_load_cobra_model():
    m = load_cobra_model(
        resource_filename("pymgpipe", "resources/models/small_sample_model.xml.gz")
    )
    assert isinstance(m, cobra.Model) and len(m.reactions) == 6091
