import cobra
import os
import pandas as pd
from pkg_resources import resource_filename
from pytest_check import check
from pymgpipe import get_abundances
from pymgpipe.modeling import build
import re


def test_build_diet_fecal():
    sample_data = [
        ["TaxaA", 0.1],
        ["TaxaB", 0.2],
        ["TaxaC", 0.3],
        ["TaxaD", 0.4]
    ]
    sample_data = pd.DataFrame(sample_data,columns=['','sample1'])
    sample_data.set_index(sample_data.columns[0],inplace=True)

    taxa_directory = resource_filename("pymgpipe", "resources/miniTaxa/")
    with check:
        assert os.path.exists(taxa_directory)

    pymgpipe_model = build(
        sample_data,
        sample='sample1',
        taxa_directory=taxa_directory,
        diet_fecal_compartments=True
    )

    with check:
        assert isinstance(pymgpipe_model, cobra.Model)
        assert (
            "fe" in pymgpipe_model.compartments and "d" in pymgpipe_model.compartments
        )

    built_abundances = get_abundances(pymgpipe_model).to_dict()["sample1"]
    true_abundances = sample_data['sample1'].to_dict()
    assert built_abundances == true_abundances


def test_build():
    sample_data = [
        ["TaxaA", 0.1],
        ["TaxaB", 0.2],
        ["TaxaC", 0.3],
        ["TaxaD", 0.4]
    ]
    sample_data = pd.DataFrame(sample_data,columns=['','sample1'])
    sample_data.set_index(sample_data.columns[0],inplace=True)

    taxa_directory = resource_filename("pymgpipe", "resources/miniTaxa/")
    with check:
        assert os.path.exists(taxa_directory)

    pymgpipe_model = build(
        sample_data,
        sample='sample1',
        taxa_directory=taxa_directory,
        diet_fecal_compartments=False
    )

    with check:
        assert isinstance(pymgpipe_model, cobra.Model)
        assert (
            "fe" not in pymgpipe_model.compartments
            and "d" not in pymgpipe_model.compartments
        )

    built_abundances = get_abundances(pymgpipe_model).to_dict()["sample1"]
    true_abundances = sample_data['sample1'].to_dict()
    assert built_abundances == true_abundances
