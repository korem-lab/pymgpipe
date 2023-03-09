import cobra
import re
from .utils import get_reactions, get_reverse_var


def remove_coupling_constraints(com):
    if isinstance(com, cobra.Model):
        com = com.solver
    c = [k for k in com.constraints if re.match(".*_cp$", k.name)]
    if len(c) == 0:
        return
    print("Removed %s coupling constraints from model!" % len(c))
    com.remove(c)


def add_coupling_constraints(com, u_const=0.01, C_const=400):
    if isinstance(com, cobra.Model):
        com = com.solver
    if 'coupled' in com.variables: # fixing some bug from old version of models
        com.remove('coupled')

    remove_coupling_constraints(com)
    biomass_rxns = {
        b.name.split("__")[-1]: b for b in get_reactions(com, regex="^biomass.*")
    }

    consts = []
    target_reactions = [
        r
        for r in com.variables
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
            consts.append(_get_coupled_upper_constraint(com,r,abundance,C_const,u_const))
        if r.lb < 0:
            consts.append(_get_coupled_lower_constraint(com,r,abundance,C_const,u_const))

    print("\nAdding coupling constraints for %s variables..." % len(consts))
    com.add(consts)
    com.update()


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