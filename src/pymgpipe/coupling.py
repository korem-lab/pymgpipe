import re
import cobra


def remove_coupling_constraints(com):
    if isinstance(com, cobra.Model):
        com = com.solver
    c = [k for k in com.constraints if re.match(".*_cp$", k.name)]
    if len(c) == 0:
        return
    print("Removed %s coupling constraints from model!" % len(c))
    com.remove(c)


def add_coupling_constraints(com, u_const=0.01, C_const=400):
    remove_coupling_constraints(com)
    biomass_rxns = {
        r.community_id: r for r in com.reactions.query(re.compile("^biomass.*"))
    }

    consts = []
    target_reactions = [
        r
        for r in com.reactions
        if not (
            r.id.startswith("EX_")
            or r.id.startswith("Diet_EX_")
            or r.id.startswith("DUt_")
            or r.id.startswith("UFEt")
            or "biomass" in r.id
            or "community" in r.id.lower()
        )
    ]
    for r in target_reactions:
        taxon = r.community_id
        abundance = com.variables[biomass_rxns[taxon].id]

        forward = r.forward_variable
        reverse = r.reverse_variable

        if forward.ub > 0:
            consts.append(
                com.solver.interface.Constraint(
                    forward - (abundance * C_const),
                    ub=u_const,
                    name="%s_cp" % forward.name,
                )
            )
        if reverse.ub > 0:
            consts.append(
                com.solver.interface.Constraint(
                    reverse - (abundance * C_const),
                    ub=u_const,
                    name="%s_cp" % reverse.name,
                )
            )

    print("\nAdding coupling constraints for %s variables..." % len(consts))
    com.solver.add(consts)
    com.solver.update()
