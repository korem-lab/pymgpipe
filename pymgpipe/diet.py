import pandas as pd
import os
import cobra
from pkg_resources import resource_listdir, resource_filename
from .utils import (
    load_dataframe,
    get_reactions,
    get_reverse_var,
    set_reaction_bounds,
    get_reaction_bounds,
)
from .io import load_model
from .logger import logger

def get_available_diets():
    """Returns all diets that come pre-packaged with pymgpipe"""
    return [f.split('.txt')[0] for f in os.listdir(resource_filename("pymgpipe", "resources/diets/"))]

def _get_adapted_diet(
    diet, essential_metabolites=None, micronutrients=None, vaginal=False, threshold=0.8
):
    if not isinstance(diet, pd.DataFrame):
        raise Exception(
            "Diet needs to be of type pd.DataFrame, received %s" % type(diet)
        )

    adapted_diet = pd.DataFrame.copy(diet)
    adapted_diet.columns = ["lb"]

    adapted_diet.replace(to_replace=r"\[e\]", value=r"(e)", regex=True, inplace=True)
    adapted_diet.replace(
        {
            "EX_adpcbl(e)": "EX_adocbl(e)",
            "EX_glc(e)": "EX_glc_D(e)",
            "EX_sbt-d(e)": "EX_sbt_D(e)",
        },
        inplace=True,
    )

    essential_metabolites = essential_metabolites if essential_metabolites is not None else (
        [
            "EX_12dgr180(e)",
            "EX_26dap_M(e)",
            "EX_2dmmq8(e)",
            "EX_2obut(e)",
            "EX_3mop(e)",
            "EX_4abz(e)",
            "EX_4hbz(e)",
            "EX_ac(e)",
            "EX_acgam(e)",
            "EX_acmana(e)",
            "EX_acnam(e)",
            "EX_ade(e)",
            "EX_adn(e)",
            "EX_adocbl(e)",
            "EX_adpcbl(e)",
            "EX_ala_D(e)",
            "EX_ala_L(e)",
            "EX_amet(e)",
            "EX_amp(e)",
            "EX_arab_D(e)",
            "EX_arab_L(e)",
            "EX_arg_L(e)",
            "EX_asn_L(e)",
            "EX_btn(e)",
            "EX_ca2(e)",
            "EX_cbl1(e)",
            "EX_cgly(e)",
            "EX_chor(e)",
            "EX_chsterol(e)",
            "EX_cit(e)",
            "EX_cl(e)",
            "EX_cobalt2(e)",
            "EX_csn(e)",
            "EX_cu2(e)",
            "EX_cys_L(e)",
            "EX_cytd(e)",
            "EX_dad_2(e)",
            "EX_dcyt(e)",
            "EX_ddca(e)",
            "EX_dgsn(e)",
            "EX_fald(e)",
            "EX_fe2(e)",
            "EX_fe3(e)",
            "EX_fol(e)",
            "EX_for(e)",
            "EX_gal(e)",
            "EX_glc_D(e)",
            "EX_gln_L(e)",
            "EX_glu_L(e)",
            "EX_gly(e)",
            "EX_glyc(e)",
            "EX_glyc3p(e)",
            "EX_gsn(e)",
            "EX_gthox(e)",
            "EX_gthrd(e)",
            "EX_gua(e)",
            "EX_h(e)",
            "EX_h2o(e)",
            "EX_h2s(e)",
            "EX_his_L(e)",
            "EX_hxan(e)",
            "EX_ile_L(e)",
            "EX_k(e)",
            "EX_lanost(e)",
            "EX_leu_L(e)",
            "EX_lys_L(e)",
            "EX_malt(e)",
            "EX_met_L(e)",
            "EX_mg2(e)",
            "EX_mn2(e)",
            "EX_mqn7(e)",
            "EX_mqn8(e)",
            "EX_nac(e)",
            "EX_ncam(e)",
            "EX_nmn(e)",
            "EX_no2(e)",
            "EX_ocdca(e)",
            "EX_ocdcea(e)",
            "EX_orn(e)",
            "EX_phe_L(e)",
            "EX_pheme(e)",
            "EX_pi(e)",
            "EX_pnto_R(e)",
            "EX_pro_L(e)",
            "EX_ptrc(e)",
            "EX_pydx(e)",
            "EX_pydxn(e)",
            "EX_q8(e)",
            "EX_rib_D(e)",
            "EX_ribflv(e)",
            "EX_ser_L(e)",
            "EX_sheme(e)",
            "EX_so4(e)",
            "EX_spmd(e)",
            "EX_thm(e)",
            "EX_thr_L(e)",
            "EX_thymd(e)",
            "EX_trp_L(e)",
            "EX_ttdca(e)",
            "EX_tyr_L(e)",
            "EX_ura(e)",
            "EX_val_L(e)",
            "EX_xan(e)",
            "EX_xyl_D(e)",
            "EX_zn2(e)",
            "EX_sucr(e)",
        ]
        if vaginal
        else [
            "EX_12dgr180(e)",
            "EX_26dap_M(e)",
            "EX_2dmmq8(e)",
            "EX_2obut(e)",
            "EX_3mop(e)",
            "EX_4abz(e)",
            "EX_4hbz(e)",
            "EX_ac(e)",
            "EX_acgam(e)",
            "EX_acmana(e)",
            "EX_acnam(e)",
            "EX_ade(e)",
            "EX_adn(e)",
            "EX_adocbl(e)",
            "EX_ala_D(e)",
            "EX_ala_L(e)",
            "EX_amet(e)",
            "EX_amp(e)",
            "EX_arab_D(e)",
            "EX_arab_L(e)",
            "EX_arg_L(e)",
            "EX_asn_L(e)",
            "EX_btn(e)",
            "EX_ca2(e)",
            "EX_cbl1(e)",
            "EX_cgly(e)",
            "EX_chor(e)",
            "EX_chsterol(e)",
            "EX_cit(e)",
            "EX_cl(e)",
            "EX_cobalt2(e)",
            "EX_csn(e)",
            "EX_cu2(e)",
            "EX_cys_L(e)",
            "EX_cytd(e)",
            "EX_dad_2(e)",
            "EX_dcyt(e)",
            "EX_ddca(e)",
            "EX_dgsn(e)",
            "EX_fald(e)",
            "EX_fe2(e)",
            "EX_fe3(e)",
            "EX_fol(e)",
            "EX_for(e)",
            "EX_gal(e)",
            "EX_glc_D(e)",
            "EX_gln_L(e)",
            "EX_glu_L(e)",
            "EX_gly(e)",
            "EX_glyc(e)",
            "EX_glyc3p(e)",
            "EX_gsn(e)",
            "EX_gthox(e)",
            "EX_gthrd(e)",
            "EX_gua(e)",
            "EX_h(e)",
            "EX_h2o(e)",
            "EX_h2s(e)",
            "EX_his_L(e)",
            "EX_hxan(e)",
            "EX_ile_L(e)",
            "EX_k(e)",
            "EX_lanost(e)",
            "EX_leu_L(e)",
            "EX_lys_L(e)",
            "EX_malt(e)",
            "EX_met_L(e)",
            "EX_mg2(e)",
            "EX_mn2(e)",
            "EX_mqn7(e)",
            "EX_mqn8(e)",
            "EX_nac(e)",
            "EX_ncam(e)",
            "EX_nmn(e)",
            "EX_no2(e)",
            "EX_ocdca(e)",
            "EX_ocdcea(e)",
            "EX_orn(e)",
            "EX_phe_L(e)",
            "EX_pheme(e)",
            "EX_pi(e)",
            "EX_pnto_R(e)",
            "EX_pro_L(e)",
            "EX_ptrc(e)",
            "EX_pydx(e)",
            "EX_pydxn(e)",
            "EX_q8(e)",
            "EX_rib_D(e)",
            "EX_ribflv(e)",
            "EX_ser_L(e)",
            "EX_sheme(e)",
            "EX_so4(e)",
            "EX_spmd(e)",
            "EX_thm(e)",
            "EX_thr_L(e)",
            "EX_thymd(e)",
            "EX_trp_L(e)",
            "EX_ttdca(e)",
            "EX_tyr_L(e)",
            "EX_ura(e)",
            "EX_val_L(e)",
            "EX_xan(e)",
            "EX_xyl_D(e)",
            "EX_zn2(e)",
            "EX_glu_D(e)",
            "EX_melib(e)",
            "EX_chtbs(e)",
            "EX_metsox_S_L(e)",
            "EX_hdca(e)",
            "EX_gam(e)",
            "EX_indole(e)",
            "EX_glcn(e)",
            "EX_coa(e)",
            "EX_man(e)",
            "EX_fum(e)",
            "EX_succ(e)",
            "EX_no3(e)",
            "EX_ins(e)",
            "EX_uri(e)",
            "EX_drib(e)",
            "EX_pime(e)",
            "EX_lac_L(e)",
            "EX_glypro(e)",
            "EX_urea(e)",
            "EX_duri(e)",
            "EX_h2(e)",
            "EX_mal_L(e)",
            "EX_tre(e)",
            "EX_orot(e)",
            "EX_glymet(e)",
            "EX_glyleu(e)",
            "EX_pydx5p(e)",
            "EX_so3(e)",
            "EX_nh4(e)",
        ]
    )

    for m in essential_metabolites:
        if (
            m not in diet.index
        ):  # this is a bug, this should be checking adapated_diet.index to account for renamed metabolites
            adapted_diet.loc[m] = [0.1]

    if not vaginal:
        unmapped_metabs = [
            "EX_asn_L(e)",
            "EX_gln_L(e)",
            "EX_crn(e)",
            "EX_elaid(e)",
            "EX_hdcea(e)",
            "EX_dlnlcg(e)",
            "EX_adrn(e)",
            "EX_hco3(e)",
            "EX_sprm(e)",
            "EX_carn(e)",
            "EX_7thf(e)",
            "EX_Lcystin(e)",
            "EX_hista(e)",
            "EX_orn(e)",
            "EX_ptrc(e)",
            "EX_creat(e)",
            "EX_cytd(e)",
            "EX_so4(e)",
        ]
        for m in unmapped_metabs:
            if m not in diet.index:
                adapted_diet.loc[m] = [50]

        if "EX_chol(e)" not in adapted_diet.index:
            adapted_diet.loc["EX_chol(e)"] = [41.251]

    micronutrients = micronutrients if micronutrients is not None else (
        [
            "EX_adocbl(e)",
            "EX_vitd2(e)",
            "EX_vitd3(e)",
            "EX_psyl(e)",
            "EX_gum(e)",
            "EX_bglc(e)",
            "EX_phyQ(e)",
            "EX_fol(e)",
            "EX_5mthf(e)",
            "EX_q10(e)",
            "EX_retinol_9_cis(e)",
            "EX_pydxn(e)",
            "EX_pydam(e)",
            "EX_pydx(e)",
            "EX_pheme(e)",
            "EX_ribflv(e)",
            "EX_thm(e)",
            "EX_avite1(e)",
            "EX_pnto_R(e)",
            "EX_na1(e)",
            "EX_cl(e)",
            "EX_k(e)",
            "EX_pi(e)",
            "EX_zn2(e)",
            "EX_cu2(e)",
        ]
        if vaginal
        else [
            "EX_adocbl(e)",
            "EX_vitd2(e)",
            "EX_vitd3(e)",
            "EX_psyl(e)",
            "EX_gum(e)",
            "EX_bglc(e)",
            "EX_phyQ(e)",
            "EX_fol(e)",
            "EX_5mthf(e)",
            "EX_q10(e)",
            "EX_retinol_9_cis(e)",
            "EX_pydxn(e)",
            "EX_pydam(e)",
            "EX_pydx(e)",
            "EX_pheme(e)",
            "EX_ribflv(e)",
            "EX_thm(e)",
            "EX_avite1(e)",
            "EX_na1(e)",
            "EX_cl(e)",
            "EX_k(e)",
            "EX_pi(e)",
            "EX_zn2(e)",
            "EX_cu2(e)",
        ]
    )

    adapted_diet.lb = -adapted_diet.lb
    adapted_diet["ub"] = threshold * adapted_diet.lb
    adapted_diet["ub"].loc[~adapted_diet.index.isin(list(diet.index))] = 0

    adapted_diet.loc[
        (adapted_diet.index.isin(set(micronutrients)))
        & (abs(adapted_diet["lb"]) <= 0.1),
        "lb",
    ] = (
        adapted_diet.lb * 100
    )
    adapted_diet.loc[
        (adapted_diet.index.isin(["EX_fol(e)", "EX_arab_L(e)"]))
        & (abs(adapted_diet["lb"]) < 1),
        "lb",
    ] = -1

    if vaginal:
        adapted_diet.loc["EX_cytd(e)", "lb"] = -10
        adapted_diet.loc["EX_ttdca(e)", "lb"] = -1

    adapted_diet.index = adapted_diet.index.str.split("\(e\)").str[0]
    return adapted_diet


# Removes any diet
def remove_diet(model):
    """Removes existing diet from model"""
    model = load_model(model)
    print("Removing diet from model...")
    for d in get_reactions(model, regex="Diet_EX_.*"):
        set_reaction_bounds(model, d, -1000, 1000)


# Finds diet current set in model
def get_diet(model):
    """Returns existing diet from model"""
    model = load_model(model)
    print("Fetching diet from model...")

    diet = []
    for f in get_reactions(model, regex="Diet_EX_.*"):
        lower, upper = get_reaction_bounds(model, f)
        diet.append({"id": f.name, "lb": lower, "ub": upper})
    return pd.DataFrame(diet)


def add_diet_to_model(
    model,
    diet,
    force_uptake=True,
    essential_metabolites=None,
    micronutrients=None,
    vaginal=False,
    threshold=0.8,
    check=True
):
    """Add pymgpipe-adapated diet to model as defined by original mgPipe paper (see README for more details)

    Args:
        model (optlang.interface.model): LP problem
        diet (pandas.DataFrame | str): Path to diet or dataframe
        essential_metabolites (list): Custom of essential metabolites  (uses pre-defined list by default)
        micronutrients (list): Custom list of micronutrients (uses pre-defined list by default)
        vaginal (bool): Whether or not this is a vaginal diet
        threshold (float): Value between 0 and 1 that defines how strict the diet constraints are (with 1 being the least strict)
        check (bool): Check whether or not this diet is feasible (can take some time depending on size of model)
    """
    model = load_model(model)

    print("\nAttempting to add diet...")
    if isinstance(diet, str) and os.path.exists(diet):
        if diet.endswith(".csv"):
            diet_df = load_dataframe(diet)
        elif diet.endswith(".txt"):
            diet_df = pd.read_csv(
                diet,
                sep="\t",
                header=0,
                index_col=0,
            )
        else:
            raise Exception(
                "Unrecognized diet file format for %s- must be .txt or .csv!" % diet
            )

        # In case of personalized diets
        if model.name in diet_df.columns:
            diet = diet_df[model.name].to_frame()
    elif isinstance(diet, str):
        try:
            diet_df = pd.read_csv(
                resource_filename("pymgpipe", "resources/diets/%s.txt" % diet),
                sep="\t",
                header=0,
                index_col=0,
            )
            print(
                "Found %s diet containing %s metabolites in resources!"
                % (diet, len(diet_df.index))
            )
        except Exception:
            logger.warning(
                "Skipping diet! Given diet `%s` not in list of available diets. Available diets are- %s"
                % (
                    diet,
                    [
                        x.split(".")[0]
                        for x in resource_listdir("pymgpipe", "resources/diets/")
                    ],
                )
            )
            return
    elif isinstance(diet, pd.DataFrame):
        print("Using custom diet with %s metabolites!" % len(diet.index))
        diet_df = diet
    else:
        logger.warning(
            "Diet not used- please pass in valid DataFrame, local file, or resource file name!"
        )
        return

    diet_reactions = get_reactions(model, regex="Diet_EX_.*")
    for d in diet_reactions:
        set_reaction_bounds(model, d, 0, 1000)

    diet_df = diet_df[diet_df.columns[0]].to_frame()

    if essential_metabolites is not None:
        print("Using custom set of essential metabolites...")
    if micronutrients is not None:
        print("Using custom set of micronutrients...")

    d = _get_adapted_diet(diet_df, essential_metabolites, micronutrients, vaginal, threshold)

    logger.info("Adding %s diet to model..." % diet)
    added = []

    diet_reactions = {
        f.name.split("[d]")[0].split("Diet_")[-1]: f for f in diet_reactions
    }
    for ex, row in d.iterrows():
        if ex in diet_reactions:
            f = diet_reactions[ex]
            if force_uptake:
                set_reaction_bounds(model, f, row.lb, row.ub)
            else:
                set_reaction_bounds(model, f, row.lb, 0)

            added.append({"id": f.name, "lb": row.lb, "ub": row.ub})

    if check:
        print("Checking diet feasibility...\n")
        model.configuration.presolve = True
        model.configuration.lp_method = 'auto'
        model.optimize()
        if model.status == "infeasible":
            logger.warning("%s is infeasible with provided diet!" % model.name)

    if len(added) == 0:
        logger.warning("Zero metabolites from diet were found within model!")

    return pd.DataFrame(added)
