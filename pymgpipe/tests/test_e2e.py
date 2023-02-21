import optlang
import os
import pandas as pd
import cobra
from pkg_resources import resource_filename
from pymgpipe import add_coupling_constraints, compute_nmpcs, build, get_abundances
from pytest_check import check
import re


def test_full_diet_fecal_compartments():
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

    pymgpipe_model = build(
        name="A test model",
        taxonomy=sample_df,
        rel_threshold=1e-6,
        solver="gurobi",
        coupling_constraints=False,
        diet_fecal_compartments=True,
    )

    add_coupling_constraints(pymgpipe_model)
    assert len([k for k in pymgpipe_model.constraints if re.match(".*_cp$", k.name)]) > 0

    built_abundances = get_abundances(pymgpipe_model).to_dict()["A test model"]
    true_abundances = sample_df.set_index("strain")["abundance"].to_dict()
    assert built_abundances == true_abundances

    nmpc_res = compute_nmpcs(
        samples=pymgpipe_model, force=True, diet_fecal_compartments=True
    )

    assert (
        os.path.exists("nmpcs.csv")
        and os.path.exists("all_fluxes.csv")
        and os.path.exists("community_objectives.csv")
    )
    os.remove("nmpcs.csv")
    os.remove("all_fluxes.csv")
    os.remove("community_objectives.csv")

    assert len(nmpc_res.nmpc) == 20


def test_full_single_compartment():
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

    pymgpipe_model = build(
        name="A test model",
        taxonomy=sample_df,
        rel_threshold=1e-6,
        solver="gurobi",
        coupling_constraints=False,
        diet_fecal_compartments=False,
    )

    add_coupling_constraints(pymgpipe_model)
    assert len([k for k in pymgpipe_model.constraints if re.match(".*_cp$", k.name)]) > 0

    built_abundances = get_abundances(pymgpipe_model).to_dict()["A test model"]
    true_abundances = sample_df.set_index("strain")["abundance"].to_dict()
    assert built_abundances == true_abundances

    nmpc_res = compute_nmpcs(
        samples=pymgpipe_model, force=True, diet_fecal_compartments=False
    )

    assert (
        os.path.exists("nmpcs.csv")
        and os.path.exists("all_fluxes.csv")
        and os.path.exists("community_objectives.csv")
    )
    os.remove("nmpcs.csv")
    os.remove("all_fluxes.csv")
    os.remove("community_objectives.csv")

    assert len(nmpc_res.nmpc) == 20
