import logging
import cobra
from optlang.symbolics import Zero
from micom.util import (
    load_model,
    clean_ids,
    compartment_id,
    COMPARTMENT_RE,
)
import re

def build(
        taxonomy,
        name,
        rel_threshold=1e-6,
        solver='gurobi',
        add_coupling_constraints=True,
        add_fecal_diet_compartments=True,
    ):
    logging.info("building new micom model {}.".format(id))
    if not solver:
        solver = [
            s
            for s in ["cplex", "gurobi", "osqp", "glpk"]
            if s in cobra.util.solver.solvers
        ][0]
    logging.info("using the %s solver." % solver)
    if solver == "glpk":
        logging.warning(
            "No QP solver found, will use GLPK. A lot of functionality "
            "in MICOM will require a QP solver :/"
        )
   
    taxonomy = taxonomy.copy()
    if "abundance" not in taxonomy.columns:
        taxonomy["abundance"] = 1
    taxonomy.abundance /= taxonomy.abundance.sum()
    logging.info(
        "{} individuals with abundances below threshold".format(
            (taxonomy.abundance <= rel_threshold).sum()
        )
    )
    taxonomy = taxonomy[taxonomy.abundance > rel_threshold]
   
    if taxonomy.id.str.contains(r"[^A-Za-z0-9_]", regex=True).any():
        logging.warning(
            "Taxa IDs contain prohibited characters and will be reformatted."
        )
        taxonomy.id = taxonomy.id.replace(r"[^A-Za-z0-9_\s]+", "_", regex=True)

    obj = Zero
    multi_species_model = cobra.Model(name=name)
    
    taxonomy.set_index('id',inplace=True)
    index = list(taxonomy.index)
    # index = track(index, description="Building") if progress else index
    logging.info('Building sample with %s taxa...'%len(index))
    for idx in index:
        row = taxonomy.loc[idx]
        model = load_model(row.file)
        suffix = "__" + idx.replace(" ", "_").strip()
        logging.info("converting IDs for {}".format(idx))
        external = cobra.medium.find_external_compartment(model)
        external = 'u'
        logging.info(
            "Identified %s as the external compartment for %s. "
            "If that is wrong you may be in trouble..." % (external, idx)
        )
        for r in model.reactions:
            r.global_id = clean_ids(r.id).replace('(e)','[u]')
            r.id = r.global_id + suffix
            r.community_id = idx
            # avoids https://github.com/opencobra/cobrapy/issues/926
            r._compartments = None
            # SBO terms may not be maintained
            if "sbo" in r.annotation:
                del r.annotation["sbo"]
        for m in model.metabolites:
            m.global_id = clean_ids(m.id).replace('[e]','[u]')
            m.id = m.global_id + suffix
            m.compartment += suffix
            m.community_id = idx
        logging.info("adding reactions for {} to community".format(idx))
        multi_species_model.add_reactions(model.reactions)
        o = multi_species_model.solver.interface.Objective.clone(
            model.objective, model=multi_species_model.solver
        )
        obj += o.expression * row.abundance
        taxa_obj = multi_species_model.problem.Constraint(
            o.expression, name="objective_" + idx, lb=0.0
        )
        multi_species_model.add_cons_vars([taxa_obj])
        _add_exchanges(
            multi_species_model,
            model.reactions,
            add_fecal_diet_compartments
        )
        multi_species_model.solver.update()  # to avoid dangling refs due to lazy add
        
    # var to track coupling constraints
    cp_var = multi_species_model.problem.Variable(name='coupled',type='binary')
    cp_var.lb = cp_var.ub =  int(add_coupling_constraints)
    multi_species_model.add_cons_vars(cp_var)

    if add_coupling_constraints:
        _add_coupling_constraints(multi_species_model)

    l_biomass = cobra.Metabolite(id='microbeBiomass[u]',compartment='u')
    multi_species_model.add_metabolites([l_biomass])

    if add_fecal_diet_compartments:
        _add_fecal_exchange(multi_species_model,l_biomass)
    else:
        _add_lumen_exchange(multi_species_model,l_biomass)

    biomass_metabs = multi_species_model.metabolites.query(re.compile('^biomass.*'))
    biomass = cobra.Reaction(
        id="communityBiomass",
        lower_bound=1,upper_bound=1,
    )
    biomass.add_metabolites({m:-float(taxonomy.abundance[m.community_id]) for m in biomass_metabs})
    biomass.add_metabolites({l_biomass:1})
    multi_species_model.add_reaction(biomass)

    multi_species_model.objective = biomass
    return multi_species_model

def _add_exchanges(
        model,
        reactions,
        fecal_diet
    ):
    """Add exchange reactions for a new model."""
    external_compartment = 'e'
    for r in reactions:
        # Some sanity checks for whether the reaction is an exchange
        ex = external_compartment + "__" + r.community_id
        # Some AGORA models label the biomass demand as exchange
        if any(bm in r.id.lower() for bm in ["biomass", "bm"]):
            r.bounds = (0,1000)
            if 'EX_' in r.id:
                model.remove_reactions([r])
            continue
        if not cobra.medium.is_boundary_type(r, "exchange", ex):
            continue
        if not r.id.lower().startswith("ex"):
            logging.warning(
                "Reaction %s seems to be an exchange " % r.id
                + "reaction but its ID does not start with 'EX_'..."
            )

        export = len(r.reactants) == 1
        met = (r.reactants + r.products)[0]
        old_compartment = compartment_id(met)
        metab_id = re.sub(
            COMPARTMENT_RE.format(old_compartment, old_compartment),
            "",
            met.global_id,
        )
        if metab_id not in model.metabolites:
            lumen_m = met.copy()
            lumen_m.id = metab_id
            lumen_m.compartment = "u"
            lumen_m.global_id = metab_id[:-3]
            lumen_m.community_id = "lumen"
            model.add_metabolites([lumen_m])

            if fecal_diet:
                _add_fecal_exchange(model,lumen_m)
                _add_diet_exchange(model, lumen_m)
            else:
                _add_lumen_exchange(model,lumen_m)
        else:
            lumen_m = model.metabolites.get_by_id(metab_id)

        coef = 1
        r.add_metabolites({lumen_m: coef if export else -coef})
        r.id = r.id.replace('EX_','IEX_')
        r.add_metabolites({k:-v for k,v in r.metabolites.items()},combine=False)
        r.lower_bound=-1000
        r.upper_bound=1000

def _add_lumen_exchange(model,met):
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

def _add_fecal_exchange(model,met):
    f_m = met.copy()
    f_m.id = met.id.replace('[u]','[fe]')
    f_m.global_id = f_m.id[:-4]
    f_m.compartment = "fe"
    f_m.community_id = "fecal"

    model.add_metabolites([f_m])

    f_ex = cobra.Reaction(
        id="EX_" + f_m.id,
        name=f_m.id + " fecal exchange",
        lower_bound=-1000,
        upper_bound=1000,
    )
    f_ex.add_metabolites({f_m:-1})
    f_ex.global_id = f_ex.id
    f_ex.community_id = "fecal"
    
    f_tr = cobra.Reaction(
        id="UFEt_" + f_m.global_id,
        name=f_m.id + " lumen fecal transport",
        lower_bound=0,
        upper_bound=1000,
    )
    f_tr.add_metabolites({met:-1,f_m:1})
    f_tr.global_id = f_tr.id
    f_tr.community_id = "lumen/fecal"

    model.add_reactions([f_ex,f_tr])

def _add_diet_exchange(model, met):
    d_m = met.copy()
    d_m.id = met.id.replace('[u]','[d]')
    d_m.global_id = d_m.id[:-3]
    d_m.compartment = "d"
    d_m.community_id = "diet"

    model.add_metabolites([d_m])

    d_ex = cobra.Reaction(
        id="EX_" + d_m.id,
        name=d_m.id + " diet exchange",
        lower_bound=-1000,
        upper_bound=1000,
    )
    d_ex.add_metabolites({d_m:-1})
    d_ex.global_id = d_ex.id
    d_ex.community_id = "diet"
    
    d_tr = cobra.Reaction(
        id="DUt_" + d_m.global_id,
        name=d_m.id + " diet lumen transport",
        lower_bound=0,
        upper_bound=1000,
    )
    d_tr.add_metabolites({met:1,d_m:-1})
    d_tr.global_id = d_tr.id
    d_tr.community_id = "diet/lumen"

    model.add_reactions([d_ex,d_tr])    

def _add_coupling_constraints(com,u_const=0.01,C_const=400):
    biomass_rxns= {r.community_id:r for r in com.reactions.query(re.compile('^biomass.*'))}

    consts = []
    target_reactions = [r for r in com.reactions if not 
        (r.id.startswith('EX_') or r.id.startswith('DUt_') or r.id.startswith('UFEt') or 
        'biomass' in r.id or 'community' in r.id.lower())
    ]
    for r in target_reactions:
        taxon = r.community_id
        abundance = com.variables[biomass_rxns[taxon].id]
        
        forward = r.forward_variable
        reverse = r.reverse_variable

        consts.append(com.solver.interface.Constraint(forward-(abundance*C_const),ub=u_const,name='%s_cp'%forward.name))
        consts.append(com.solver.interface.Constraint(reverse-(abundance*C_const),ub=u_const,name='%s_cp'%reverse.name))
    
    com.solver.add(consts)
    com.solver.update()