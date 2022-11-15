import re
import logging
import cobra

def remove_coupling_constraints(com):
    if isinstance(com,cobra.Model):
        com=com.solver
    c = [k for k in com.constraints if re.match(".*_cp$",k.name)]
    if len(c) == 0:
        logging.info('No coupling constraints to remove!')
        return
    com.remove(c)

def add_coupling_constraints(com,u_const=0.01,C_const=400):
    if isinstance(com,cobra.Model):
        com=com.solver
    remove_coupling_constraints(com)
    biomass_rxns = {r.name.split('_pan')[1]:r for r in com.variables if re.match('.*biomass.*$',r.name) and 'reverse' not in r.name}

    consts = []
    target_reactions = [k for k in com.variables if re.match('^.*pan((?!biomass).)*$',k.name)]
    for v in target_reactions:
        taxon=v.name.split('_pan')[1].split('_reverse')[0]
        abundance = biomass_rxns[taxon]

        if v.ub > 0:
            consts.append(com.interface.Constraint(v-(abundance*C_const),ub=u_const,name='%s_cp'%v.name))
    
    com.add(consts)
    com.update()