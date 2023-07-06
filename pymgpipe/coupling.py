import cobra
import re
from .utils import get_reactions, get_reverse_var


def remove_coupling_constraints(model):
    """Removing coupling constraints to community-level models

    Args:
        model (optlang.interface.Model): LP problem

    Notes:
        Removes all coupling constraints from model (if they exist)
    """
    if isinstance(model, cobra.Model):
        model = model.solver
    c = [k for k in model.constraints if re.match(".*_cp$", k.name)]
    if len(c) == 0:
        return
    print("Removed %s coupling constraints from model!" % len(c))
    model.remove(c)


def add_coupling_constraints(model, u_const=0.01, C_const=400):
    """Adding coupling constraints to community-level models

    Args:
        model (optlang.interface.Model): LP problem
        u_const (float): U flexibility constant
        C_const (float): C flexibility constant

    Notes:
        Coupling constraints are essentially a set of lower/upper bound constraints that scale reaction fluxes to the relative abundance of the microbe they belong to, like so-
        Imagine a reaction `-1000 <= My_reaction_taxa_A <= 1000`. Adding coupling constraints would redefine the bounds like so-

        `-C_const - (u_counts * <abundance of taxa A>) <= My_reaction_taxa_A <= C_const + (u_counts * <abundance of taxa A>)`
    """
    if isinstance(model, cobra.Model):
        model = model.solver
    if 'coupled' in model.variables: # fixing some bug from old version of models
        model.remove('coupled')

    remove_coupling_constraints(model)
    biomass_rxns = {
        b.name.split("__")[-1]: b for b in get_reactions(model, regex="^biomass.*")
    }

    consts = []
    target_reactions = [
        r
        for r in model.variables
        if not (
            r.name.startswith("EX_")
            or r.name.startswith("Diet_EX_")
            or r.name.startswith("DUt_")
            or r.name.startswith("UFEt")
            or "biomass" in r.name.lower()
            or "community" in r.name.lower()
        )
    ]
    for r in target_reactions:
        taxon = r.name.split('_reverse')[0].split("__")[-1]
        if taxon.startswith('_'):
            taxon = taxon[1:]
        try:
            abundance = biomass_rxns[taxon]
        except:
            raise Exception('Issue parsing taxon from reaction %s'%r.name)

        if r.ub > 0:
            consts.append(_get_coupled_upper_constraint(model,r,abundance,C_const,u_const))
        if r.lb < 0:
            consts.append(_get_coupled_lower_constraint(model,r,abundance,C_const,u_const))

    print("\nAdding coupling constraints for %s variables..." % len(consts))
    model.add(consts)
    model.update()


def _get_coupled_upper_constraint(model, v, abundance, C_const, u_const):
    return model.interface.Constraint(
        v - (abundance * C_const),
        ub=u_const,
        name="%s_cp" % v.name,
    )

def _get_coupled_lower_constraint(model, v, abundance, C_const, u_const):
    return model.interface.Constraint(
        v + (abundance * C_const),
        lb=-u_const,
        name="%s_l_cp" % v.name,
    )