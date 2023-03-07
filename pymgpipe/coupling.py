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
            or "reverse" in r.name.lower()
        )
    ]
    for r in target_reactions:
        taxon = r.name.split("__")[-1]
        abundance = biomass_rxns[taxon]

        forward = r
        forward_const = com.interface.Constraint(
            forward - (abundance * C_const),
            ub=u_const,
            name="%s_cp" % forward.name,
        )
        if forward.ub > 0:
            consts.append(forward_const)

        try:
            reverse = get_reverse_var(com, forward)
            reverse_const = com.interface.Constraint(
                reverse - (abundance * C_const),
                ub=u_const,
                name="%s_cp" % reverse.name,
            )
            if reverse.ub > 0:
                consts.append(reverse_const)
        except:
            # Catch if no reverse variables in model
            reverse_const = com.interface.Constraint(
                forward + (abundance * C_const),
                lb=-u_const,
                name="%s_l_cp" % forward.name,
            )
            if forward.lb < 0:
                consts.append(reverse_const)
                
    print("\nAdding coupling constraints for %s variables..." % len(consts))
    com.add(consts)
    com.update()
