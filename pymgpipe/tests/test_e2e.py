import optlang
import os
import pandas as pd
import numpy as np
import cobra
from pkg_resources import resource_filename
from pymgpipe import (
    add_coupling_constraints,
    compute_nmpcs,
    get_abundances,
    build_models,
    remove_reverse_vars,
    load_dataframe,
    get_reverse_id,
)
from pymgpipe.build import _build
from pytest_check import check
import re
import tempfile


def test_build_models():
    samples = ["sample%s" % i for i in range(5)]
    cov = pd.DataFrame(columns=samples, index=["TaxaA", "TaxaB", "TaxaC", "TaxaD"])

    for t in cov.columns:
        cov[t] = np.random.dirichlet(np.ones(4), size=1)[0]

    with tempfile.TemporaryDirectory() as tmpdirname:
        build_models(
            coverage_file=cov,
            taxa_dir=resource_filename("pymgpipe", "resources/miniTaxa/"),
            threads=-1,
            parallel=True,
            out_dir=tmpdirname,
            diet_fecal_compartments=True,
            coupling_constraints=True,
            compute_metrics=True,
            compress=True,
        )
        models_out = os.listdir(tmpdirname + "/models/")
        problems_out = os.listdir(tmpdirname + "/problems/")

        assert len(models_out) == 5 and len(problems_out) == 5

        assert (
            os.path.exists(tmpdirname + "/metabolic_diversity.png")
            and os.path.exists(tmpdirname + "/reaction_abundance.csv")
            and os.path.exists(tmpdirname + "/reaction_content.csv")
            and os.path.exists(tmpdirname + "/sample_label_conversion.csv")
        )


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

    pymgpipe_model = _build(
        name="A test model",
        taxonomy=sample_df,
        rel_threshold=1e-6,
        solver="gurobi",
        diet_fecal_compartments=True,
    )

    add_coupling_constraints(pymgpipe_model)
    assert (
        len([k for k in pymgpipe_model.constraints if re.match(".*_cp$", k.name)]) > 0
    )

    built_abundances = get_abundances(pymgpipe_model).to_dict()["A test model"]
    true_abundances = sample_df.set_index("strain")["abundance"].to_dict()
    assert built_abundances == true_abundances

    nmpc_res = compute_nmpcs(
        samples=pymgpipe_model, force=True, diet_fecal_compartments=True, objective_percent=None
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

    pymgpipe_model = _build(
        name="A test model",
        taxonomy=sample_df,
        rel_threshold=1e-6,
        solver="gurobi",
        diet_fecal_compartments=False,
    )

    add_coupling_constraints(pymgpipe_model)
    assert (
        len([k for k in pymgpipe_model.constraints if re.match(".*_cp$", k.name)]) > 0
    )

    built_abundances = get_abundances(pymgpipe_model).to_dict()["A test model"]
    true_abundances = sample_df.set_index("strain")["abundance"].to_dict()
    assert built_abundances == true_abundances

    nmpc_res = compute_nmpcs(
        samples=pymgpipe_model, force=True, diet_fecal_compartments=False, objective_percent=None
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


def test_remove_variables():
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

    pymgpipe_model = _build(
        name="A test model",
        taxonomy=sample_df,
        rel_threshold=1e-6,
        solver="gurobi",
        diet_fecal_compartments=True,
    )
    some_var = pymgpipe_model.variables[100]
    reverse_var_id = get_reverse_id(some_var.name)

    res1 = compute_nmpcs(samples=pymgpipe_model, write_to_file=False, threads=-1, objective_percent=None).nmpc
    remove_reverse_vars(pymgpipe_model, hard_remove=False)
    
    assert pymgpipe_model.variables[reverse_var_id].lb == 0 and pymgpipe_model.variables[reverse_var_id].ub == 0
    
    res2 = compute_nmpcs(samples=pymgpipe_model, write_to_file=False, threads=-1, objective_percent=None).nmpc
    remove_reverse_vars(pymgpipe_model, hard_remove=True)
    
    assert reverse_var_id not in pymgpipe_model.variables
    
    res3 = compute_nmpcs(samples=pymgpipe_model, write_to_file=False, threads=-1, objective_percent=None).nmpc
    
    assert (
        len(_compare(res1, res2)) == 0
        and len(_compare(res1, res3)) == 0
    )


def _compare(first, second, threshold=1e-10):
    first = load_dataframe(first)
    second = load_dataframe(second)
    bad = []
    for i, row in first.iterrows():
        second_row = second.loc[i].to_dict()
        for x in row.to_dict():
            first_x = row[x]
            second_x = second_row[x]
            if abs(first_x - second_x) > threshold:
                bad.append((i, x, abs(first_x - second_x)))
    return bad
