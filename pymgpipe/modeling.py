import cobra
import re
import os
import time
from optlang.symbolics import Zero
from cobra.medium import is_boundary_type
from .io import load_cobra_model, UnsupportedSolverException
from .utils import load_dataframe
from .logger import logger

def build(
    abundances,
    sample,
    taxa_directory,
    threshold=1e-6,
    diet_fecal_compartments=True,
    solver="gurobi",
):
    """Build community COBRA model using mgpipe-like compartments and constraints.

    This function is pymgpipe's main model building function, can can be used to build a single model (with no further modifications)

    Args:
        abundances (pandas.DataFrame | str): Abundance matrix with taxa as rows and samples as columns
        sample (str): Label corresponding to the sample you want to build (needs to match up to column name in abundance matrix)
        taxa_directory (str): Directory containing individual strain/species taxa models (file names corresponding to index of coverage matrix)
        threshold (float): Abundance threshold, any taxa with an abundance less than this value will be left out and abundances will be re-normalized
        diet_fecal_compartments (bool): Build models with mgpipe's diet/fecal compartmentalization, defaults to False
        solver (str): LP solver (gurobi or cplex) used to solve models, defaults to gurobi

    """
    abundances = load_dataframe(abundances)
    assert sample in abundances.columns, 'Sample %s not found in abundance matrix!'%sample 

    if not os.path.exists(taxa_directory):
        raise Exception('Taxa directory %s not found!'%taxa_directory)

    if solver not in ['gurobi','cplex']:
        raise UnsupportedSolverException

    sample_abundances = abundances[sample]
    sample_abundances = sample_abundances[sample_abundances != 0]

    if threshold is not None:
        sample_abundances = sample_abundances[sample_abundances > threshold]

    sample_abundances = sample_abundances / sample_abundances.sum()

    existing_taxa_files = {
        t.split("/")[-1].split(".")[0]: os.path.join(taxa_directory,t) for t in os.listdir(taxa_directory)
    }
    missing = [t for t in sample_abundances.index if t not in existing_taxa_files]
    if len(missing) > 0:
        logger.warning('Could not find associated models for %s taxa- %s\nRemoving missing taxa and renormalizing abundances.'%(len(missing),missing))
    
        sample_abundances.drop(missing, inplace=True)
        sample_abundances = sample_abundances / sample_abundances.sum()

    print('Building community model for %s with %s unique taxa...\n'%(sample,len(sample_abundances.index)))
    start = time.time()
    community_model = cobra.Model(name=sample)
    community_model.solver = solver 

    for taxon in sample_abundances.index:
        model = load_cobra_model(existing_taxa_files[taxon])
        taxon = taxon.replace(' ','_')
        ex_metabolites = [m for m in model.metabolites if '[e]' in m.id]
        missing = []
        for ex in ex_metabolites:
            if 'EX_%s(e)'%ex.id.split('[e]')[0] not in model.reactions:                
                missing.append(_get_missing_exchange(ex))
        if len(missing) > 0:
            print('Adding %s missing exchange reaction(s) to %s!'%(len(missing),taxon))
            model.add_reactions(missing)

        # -- Reactions --
        for r in list(model.reactions):
            r.id = _remove_non_alphanumeric(r.id+'__'+taxon).replace("(e)", "[u]")

        # -- Metabolites -- 
        for m in list(model.metabolites):
            m.id = _remove_non_alphanumeric(m.id+'__'+taxon).replace("[e]", "[u]")
            m.compartment = _remove_non_alphanumeric(m.compartment+'__'+taxon).replace('e__','u__')
        
        community_model.add_reactions(model.reactions)
        community_model.solver.update()

    print('Adding exchange reactions...\n')
    _add_exchanges(community_model, diet_fecal_compartments)

    microbe_biomass = cobra.Metabolite(id="microbeBiomass[u]", compartment="u")
    community_model.add_metabolites([microbe_biomass])

    if diet_fecal_compartments:
        _add_fecal_exchange(community_model, microbe_biomass)
    else:
        _add_lumen_exchange(community_model, microbe_biomass)

    biomass_metabs = community_model.metabolites.query(re.compile("^biomass.*", re.IGNORECASE))
    biomass = cobra.Reaction(
        id="communityBiomass",
        lower_bound=0.4,
        upper_bound=1,
    )
    biomass.add_metabolites(
        {m: -float(sample_abundances[m.id.split('__')[-1]]) for m in biomass_metabs}
    )
    biomass.add_metabolites({microbe_biomass: 1})
    community_model.add_reactions([biomass])

    print('Setting community objective to biomass reaction...')
    print(biomass.id+': '+biomass.reaction)

    community_model.objective = biomass
    community_model.objective_direction = 'max'
    community_model.solver.update()

    elapsed = time.time() - start
    print('\n-----------------------------------')
    logger.info('Finished building %s in %.2f minutes!'%(sample,elapsed/60))
    return community_model

def _add_exchanges(model, diet_fecal_compartments):
    for r in list(model.reactions):
        r_taxon = r.id.split('__')[-1]

        # -- Handling non-exchanges --
        if r.id.startswith("DM_"):
            r.lower_bound = 0
        if r.id.startswith("sink_"):
            r.lower_bound = -1
        if 'biomass' in r.id.lower():
            if 'EX_' in r.id:
                model.remove_reactions(r.id)
            r.bounds = (0, 1000)
        
        # -- Handling exchanges --
        if not is_boundary_type(r, 'exchange', 'u__'+r_taxon):
            continue 

        metab = (r.reactants + r.products)[0]
        lumen_id = metab.id.split('__'+r_taxon)[0]
        if lumen_id not in model.metabolites:
            lumen_metab = metab.copy()
            lumen_metab.id = lumen_id
            lumen_metab.compartment = 'u'
            
            model.add_metabolites([lumen_metab])

            if diet_fecal_compartments:
                _add_fecal_exchange(model, lumen_metab)
                _add_diet_exchange(model, lumen_metab)
            else:
                _add_lumen_exchange(model, lumen_metab)
        else:
            lumen_metab = model.metabolites.get_by_id(lumen_id)
        
        r.add_metabolites({lumen_metab: 1 if len(r.reactants) == 1 else -1})
        r.id = r.id.replace("EX_", "IEX_")
        r.add_metabolites({k: -v for k, v in r.metabolites.items()}, combine=False)
        r.lower_bound = -1000
        r.upper_bound = 1000
    model.solver.update()


def _add_lumen_exchange(model, met):
    ex_medium = cobra.Reaction(
        id="EX_" + met.id,
        name=met.id + " lumen exchange",
        lower_bound=-1000,
        upper_bound=1000,
    )
    ex_medium.add_metabolites({met: -1})
    ex_medium.global_id = ex_medium.id
    ex_medium.community_id = "lumen"
    model.add_reactions([ex_medium])


def _add_fecal_exchange(model, met):
    f_m = met.copy()
    f_m.id = met.id.replace("[u]", "[fe]")
    f_m.global_id = f_m.id[:-4]
    f_m.compartment = "fe"
    f_m.community_id = "fecal"

    model.add_metabolites([f_m])

    f_ex = cobra.Reaction(
        id="EX_" + f_m.id,
        name=f_m.id + " fecal exchange",
        lower_bound=-1000,
        upper_bound=1000000,
    )
    f_ex.add_metabolites({f_m: -1})
    f_ex.global_id = f_ex.id
    f_ex.community_id = "fecal"

    f_tr = cobra.Reaction(
        id="UFEt_" + f_m.global_id,
        name=f_m.id + " lumen fecal transport",
        lower_bound=0,
        upper_bound=1000000,
    )
    f_tr.add_metabolites({met: -1, f_m: 1})
    f_tr.global_id = f_tr.id
    f_tr.community_id = "lumen/fecal"

    model.add_reactions([f_ex, f_tr])


def _add_diet_exchange(model, met):
    d_m = met.copy()
    d_m.id = met.id.replace("[u]", "[d]")
    d_m.global_id = d_m.id[:-3]
    d_m.compartment = "d"
    d_m.community_id = "diet"

    model.add_metabolites([d_m])

    d_ex = cobra.Reaction(
        id="Diet_EX_" + d_m.id,
        name=d_m.id + " diet exchange",
        lower_bound=-1000,
        upper_bound=1000,
    )
    d_ex.add_metabolites({d_m: -1})
    d_ex.global_id = d_ex.id
    d_ex.community_id = "diet"

    d_tr = cobra.Reaction(
        id="DUt_" + d_m.global_id,
        name=d_m.id + " diet lumen transport",
        lower_bound=0,
        upper_bound=1000000,
    )
    d_tr.add_metabolites({met: 1, d_m: -1})
    d_tr.global_id = d_tr.id
    d_tr.community_id = "diet/lumen"

    model.add_reactions([d_ex, d_tr])

def _get_missing_exchange(metab):
    ex = cobra.Reaction(
        id='EX_%s(e)'%metab.id.split('[e]')[0],
        name='%s exchange'%metab.name,
        lower_bound=-1000,
        upper_bound=1000,
        subsystem='Exchange/demand reaction')
    ex.add_metabolites({metab:-1})
    return ex

def _remove_non_alphanumeric(string):
    return re.sub(r'[^\w[\]()]+', '', string)