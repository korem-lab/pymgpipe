import optlang
import os
import pandas as pd
from pymgpipe import get_reactions, regularFVA, compute_nmpcs


def test_regularFVA(mini_optlang_model):
    ex_reactions = [r.name for r in get_reactions(mini_optlang_model, regex="EX_.*")]
    fva_res = regularFVA(mini_optlang_model, ex_reactions)

    assert set(ex_reactions) == set(fva_res.index.to_list())


def test_nmpc(mini_optlang_model):
    nmpc_res = compute_nmpcs(samples=mini_optlang_model, force=True)

    assert (
        os.path.exists("nmpcs.csv")
        and os.path.exists("all_fluxes.csv")
        and os.path.exists("community_objectives.csv")
    )
    os.remove("nmpcs.csv")
    os.remove("all_fluxes.csv")
    os.remove("community_objectives.csv")

    assert len(nmpc_res.nmpc) == 20


def test_nmpc_multisample(mini_optlang_model):
    second_sample = optlang.Model.clone(mini_optlang_model)
    second_sample.name = "A second sample"

    nmpc_res = compute_nmpcs(samples=[mini_optlang_model, second_sample], force=True)

    assert (
        os.path.exists("nmpcs.csv")
        and os.path.exists("all_fluxes.csv")
        and os.path.exists("community_objectives.csv")
    )
    os.remove("nmpcs.csv")
    os.remove("all_fluxes.csv")
    os.remove("community_objectives.csv")

    assert (
        nmpc_res.nmpc["mini_model"]
        .round(6)
        .equals(nmpc_res.nmpc["A second sample"].round(6))
    )
