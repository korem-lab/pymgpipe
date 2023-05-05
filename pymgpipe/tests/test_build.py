import cobra
import os
import pandas as pd
from pkg_resources import resource_filename
from pytest_check import check
from pymgpipe import get_abundances
from pymgpipe.build import build
import re


def test_build_diet_fecal():
    sample_data = [
        ["mc1", 0.1, "TaxaA"],
        ["mc1", 0.2, "TaxaB"],
        ["mc1", 0.3, "TaxaC"],
        ["mc1", 0.4, "TaxaD"],
    ]

    sample_df = pd.DataFrame(sample_data, columns=["sample_id", "abundance", "strain"])
    sample_df["id"] = sample_df["strain"]
    sample_df["file"] = (
        resource_filename("pymgpipe", "resources/miniTaxa/") + sample_df.id + ".xml.gz"
    )

    with check:
        assert os.path.exists(sample_df.file[0])

    pymgpipe_model = build(
        name="A test model",
        taxonomy=sample_df,
        rel_threshold=1e-6,
        solver="gurobi",
        diet_fecal_compartments=True,
    )

    with check:
        assert isinstance(pymgpipe_model, cobra.Model)
        assert (
            "fe" in pymgpipe_model.compartments and "d" in pymgpipe_model.compartments
        )

    built_abundances = get_abundances(pymgpipe_model).to_dict()["A test model"]
    true_abundances = sample_df.set_index("strain")["abundance"].to_dict()
    assert built_abundances == true_abundances


def test_build():
    sample_data = [
        ["mc1", 0.1, "TaxaA"],
        ["mc1", 0.2, "TaxaB"],
        ["mc1", 0.3, "TaxaC"],
        ["mc1", 0.4, "TaxaD"],
    ]

    sample_df = pd.DataFrame(sample_data, columns=["sample_id", "abundance", "strain"])
    sample_df["id"] = sample_df["strain"]
    sample_df["file"] = (
        resource_filename("pymgpipe", "resources/miniTaxa/") + sample_df.id + ".xml.gz"
    )

    with check:
        assert os.path.exists(sample_df.file[0])

    pymgpipe_model = build(
        name="A test model",
        taxonomy=sample_df,
        rel_threshold=1e-6,
        solver="gurobi",
        diet_fecal_compartments=False,
    )

    with check:
        assert isinstance(pymgpipe_model, cobra.Model)
        assert (
            "fe" not in pymgpipe_model.compartments
            and "d" not in pymgpipe_model.compartments
        )

    built_abundances = get_abundances(pymgpipe_model).to_dict()["A test model"]
    true_abundances = sample_df.set_index("strain")["abundance"].to_dict()
    assert built_abundances == true_abundances
