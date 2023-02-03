import pytest
import optlang
import os
from pymgpipe import get_reactions, regularFVA, compute_nmpcs, load_model


def test_fva(small_optlang_model):
    ex_reactions = [
        r.name for r in get_reactions(small_optlang_model, regex="EX_.*")[:20]
    ]
    fva_res = regularFVA(small_optlang_model, ex_reactions)

    assert fva_res.index.to_list() == ex_reactions


def test_nmpc(small_optlang_model):
    nmpc_res = compute_nmpcs(samples=small_optlang_model, force=True)

    assert (
        os.path.exists("nmpcs.csv")
        and os.path.exists("all_fluxes.csv")
        and os.path.exists("community_objectives.csv")
    )
    os.remove("nmpcs.csv")
    os.remove("all_fluxes.csv")
    os.remove("community_objectives.csv")

    assert len(nmpc_res.nmpc) == 215


def test_nmpc_multisample(small_optlang_model):
    second_sample = optlang.Model.clone(small_optlang_model)
    second_sample.name = "A second sample"

    nmpc_res = compute_nmpcs(samples=[small_optlang_model, second_sample], force=True)

    assert (
        os.path.exists("nmpcs.csv")
        and os.path.exists("all_fluxes.csv")
        and os.path.exists("community_objectives.csv")
    )
    os.remove("nmpcs.csv")
    os.remove("all_fluxes.csv")
    os.remove("community_objectives.csv")

    assert (
        nmpc_res.nmpc["small_sample_model"]
        .round(6)
        .equals(nmpc_res.nmpc["A second sample"].round(6))
    )
