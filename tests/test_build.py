from pytest_check import check
from pymgpipe import *
import os 

def test_build():
    sample_data = [
        ['mc1',0.01534,'panLactobacillus_iners'],
        ['mc1',0.25734,'panGardnerella_vaginalis'],
        ['mc1',0.26341,'panLactobacillus_crispatus'],
        ['mc1',0.003123,'panLactobacillus_jensenii'],
        ['mc1',0.15342,'panMegasphaera_elsdenii'],
        ['mc1',0.307367,'panAtopobium_vaginae']
    ]
    sample_df = pd.DataFrame(sample_data, columns=['sample_id','abundance','strain'])
    sample_df['id']=sample_df['strain']
    sample_df['file']= resource_filename('pymgpipe','resources/taxaModels/') + sample_df.id + '.xml.gz'

    with check:
        assert(os.path.exists(sample_df.file[0]))

    pymgpipe_model = build(
        name='A test model',
        taxonomy=sample_df,
        rel_threshold=1e-6,
        solver='gurobi',
        coupling_constraints=True,
        diet_fecal_compartments=True
    )

    with check:
        assert(isinstance(pymgpipe_model, cobra.Model))

    built_abundances = get_abundances(pymgpipe_model).to_dict()['A test model']
    true_abundances = sample_df.set_index('strain')['abundance'].to_dict()
    assert built_abundances == true_abundances