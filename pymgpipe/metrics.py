from .utils import get_abundances
import cobra

def _compute_diversity_metrics(model):
    assert isinstance(model, cobra.Model), '`model` needs to be COBRA model to compute diversity metrics.'
    print('Computing metrics for %s...'%model.name)

    taxa = get_abundances(model)[model.name].to_dict()
    unique_reactions = set(r.id.split('__')[0] for r in model.reactions if 
        'EX' not in r.id and 
        'biomass' not in r.id and
        'UFEt_' not in r.id and 
        'DUt_' not in r.id and
        'community' not in r.id and
        'sink' not in r.id
    )
    rxn_abundance = {r:0 for r in unique_reactions}
    for r in model.reactions:
        try:
            rxn_id = r.id.split('__')[0]
            if not rxn_id in unique_reactions:
                continue 

            rxn_taxa = r.id.split('__')[1]
            rxn_taxa = rxn_taxa.split('_')[1] if rxn_taxa.startswith('_') else rxn_taxa # small fix for weird naming bug
            rxn_abundance[rxn_id] += taxa[rxn_taxa]
        except:
            pass

    to_return = {}
    to_return['sample']=model.name
    to_return['taxa']= list(taxa.keys())
    to_return['unique_reactions']=unique_reactions
    to_return['reaction_abundance'] = rxn_abundance
    return to_return

