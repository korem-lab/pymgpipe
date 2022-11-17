import pandas as pd
from pkg_resources import resource_listdir,resource_filename
from .utils import *
from .utils import Constants
import logging

def get_adapted_diet(diet,essential_metabolites=None):
    if not isinstance(diet,pd.DataFrame):
        raise Exception('Diet needs to be of type pd.DataFrame, received %s'%type(diet))
  
    diet.columns = ['metab','lb']
    diet.set_index(['metab'],inplace=True)
    adapted_diet = pd.DataFrame.copy(diet)

    adapted_diet.replace(to_replace=r'\[e\]',value=r'(e)',regex=True,inplace=True)
    adapted_diet.replace({
        'EX_adpcbl(e)':'EX_adocbl(e)',
        'EX_glc(e)':'EX_glc_D(e)',
        'EX_sbt-d(e)':'EX_sbt_D(e)',
    },inplace=True)

    essential_metabolites = essential_metabolites if essential_metabolites is not None else ['EX_12dgr180(e)', 'EX_26dap_M(e)', 'EX_2dmmq8(e)', 'EX_2obut(e)', 'EX_3mop(e)', 'EX_4abz(e)', 'EX_4hbz(e)', 'EX_ac(e)', 'EX_acgam(e)', 'EX_acmana(e)', 'EX_acnam(e)', 'EX_ade(e)', 'EX_adn(e)', 'EX_adocbl(e)', 'EX_adpcbl(e)', 'EX_ala_D(e)', 'EX_ala_L(e)', 'EX_amet(e)', 'EX_amp(e)', 'EX_arab_D(e)', 'EX_arab_L(e)', 'EX_arg_L(e)', 'EX_asn_L(e)', 'EX_btn(e)', 'EX_ca2(e)', 'EX_cbl1(e)', 'EX_cgly(e)', 'EX_chor(e)', 'EX_chsterol(e)', 'EX_cit(e)', 'EX_cl(e)', 'EX_cobalt2(e)', 'EX_csn(e)', 'EX_cu2(e)', 'EX_cys_L(e)', 'EX_cytd(e)', 'EX_dad_2(e)', 'EX_dcyt(e)', 'EX_ddca(e)', 'EX_dgsn(e)', 'EX_fald(e)', 'EX_fe2(e)', 'EX_fe3(e)', 'EX_fol(e)', 'EX_for(e)', 'EX_gal(e)', 'EX_glc_D(e)', 'EX_gln_L(e)', 'EX_glu_L(e)', 'EX_gly(e)', 'EX_glyc(e)', 'EX_glyc3p(e)', 'EX_gsn(e)', 'EX_gthox(e)', 'EX_gthrd(e)', 'EX_gua(e)', 'EX_h(e)', 'EX_h2o(e)', 'EX_h2s(e)', 'EX_his_L(e)', 'EX_hxan(e)', 'EX_ile_L(e)', 'EX_k(e)', 'EX_lanost(e)', 'EX_leu_L(e)', 'EX_lys_L(e)', 'EX_malt(e)', 'EX_met_L(e)', 'EX_mg2(e)', 'EX_mn2(e)', 'EX_mqn7(e)', 'EX_mqn8(e)', 'EX_nac(e)', 'EX_ncam(e)', 'EX_nmn(e)', 'EX_no2(e)', 'EX_ocdca(e)', 'EX_ocdcea(e)', 'EX_orn(e)', 'EX_phe_L(e)', 'EX_pheme(e)', 'EX_pi(e)', 'EX_pnto_R(e)', 'EX_pro_L(e)', 'EX_ptrc(e)', 'EX_pydx(e)', 'EX_pydxn(e)', 'EX_q8(e)', 'EX_rib_D(e)', 'EX_ribflv(e)', 'EX_ser_L(e)', 'EX_sheme(e)', 'EX_so4(e)', 'EX_spmd(e)', 'EX_thm(e)', 'EX_thr_L(e)', 'EX_thymd(e)', 'EX_trp_L(e)', 'EX_ttdca(e)', 'EX_tyr_L(e)', 'EX_ura(e)', 'EX_val_L(e)', 'EX_xan(e)', 'EX_xyl_D(e)', 'EX_zn2(e)', 'EX_glu_D(e)', 'EX_melib(e)', 'EX_chtbs(e)', 'EX_metsox_S_L(e)', 'EX_hdca(e)', 'EX_gam(e)', 'EX_indole(e)', 'EX_glcn(e)', 'EX_coa(e)', 'EX_man(e)', 'EX_fum(e)', 'EX_succ(e)', 'EX_no3(e)', 'EX_ins(e)', 'EX_uri(e)', 'EX_drib(e)', 'EX_pime(e)', 'EX_lac_L(e)', 'EX_glypro(e)', 'EX_urea(e)', 'EX_duri(e)', 'EX_h2(e)', 'EX_mal_L(e)', 'EX_tre(e)', 'EX_orot(e)', 'EX_glymet(e)', 'EX_glyleu(e)', 'EX_sucr(e)', 'EX_sbt_D(e)', 'EX_akg(e)', 'EX_glytyr(e)', 'EX_glypro(e)', 'EX_glyphe(e)', 'EX_glymet(e)', 'EX_glyleu(e)', 'EX_glyglu(e)', 'EX_glygln(e)', 'EX_glyasp(e)', 'EX_glyasn(e)', 'EX_alagly(e)', 'EX_gam(e)', 'EX_mantr(e)']
    for m in essential_metabolites:
        if m not in diet.index: # this is a bug, this should be checking adapated_diet.index to account for renamed metabolites
            adapted_diet.loc[m]=[0.1]

    unmapped_metabs = ['EX_asn_L(e)','EX_gln_L(e)','EX_crn(e)','EX_elaid(e)','EX_hdcea(e)','EX_dlnlcg(e)','EX_adrn(e)','EX_hco3(e)','EX_sprm(e)','EX_carn(e)','EX_7thf(e)','EX_Lcystin(e)','EX_hista(e)','EX_orn(e)','EX_ptrc(e)','EX_creat(e)','EX_cytd(e)','EX_so4(e)'];
    for m in unmapped_metabs:
        if m not in diet.index:
            adapted_diet.loc[m]=[50]

    if not 'EX_chol(e)' in adapted_diet.index:
        adapted_diet.loc['EX_chol(e)']=[41.251]
    
    micro_nutrients = ['EX_adocbl(e)','EX_vitd2(e)','EX_vitd3(e)','EX_psyl(e)','EX_gum(e)','EX_bglc(e)','EX_phyQ(e)','EX_fol(e)','EX_5mthf(e)','EX_q10(e)','EX_retinol_9_cis(e)','EX_pydxn(e)','EX_pydam(e)','EX_pydx(e)','EX_pheme(e)','EX_ribflv(e)','EX_thm(e)','EX_avite1(e)','EX_na1(e)','EX_cl(e)','EX_k(e)','EX_pi(e)','EX_zn2(e)','EX_cu2(e)']

    adapted_diet.lb = -adapted_diet.lb
    adapted_diet['ub'] = 0.8 * adapted_diet.lb
    adapted_diet['ub'].loc[~adapted_diet.index.isin(list(diet.index))] = 0

    adapted_diet.loc[(adapted_diet.index.isin(set(micro_nutrients))) & (abs(adapted_diet['lb']) <= 0.1),'lb']=adapted_diet.lb*100
    adapted_diet.loc[(adapted_diet.index.isin(['EX_fol(e)','EX_arab_L(e)'])) & (abs(adapted_diet['lb']) < 1),'lb']=-1

    adapted_diet.index = adapted_diet.index.str.split('\(e\)').str[0]
    return adapted_diet

def remove_diet(model):
    print('Removing diet from model...')
    for m in model.reactions:
        if m.id.startswith('Diet_EX_'):
            m.bounds = (-1000,1000)

def add_diet_to_model(
    model,
    diet
):
    print('\nAttempting to add diet...')
    if isinstance(diet,str) and os.path.exists(diet):
        diet_df = pd.read_csv(diet,index_col=0)
    elif isinstance(diet,str):
        try:
            diet_df = pd.read_csv(resource_filename('pymgpipe','resources/diets/%s.txt'%diet), sep="\t", header=0)
            print('Found %s diet containing %s metabolites in resources!'%(diet,len(diet_df.index))) 
        except:
            logging.warn('Skipping diet! Given diet `%s` not in list of available diets. Available diets are- %s'%(diet,[x.split('.')[0] for x in resource_listdir('pymgpipe','resources/diets/')]))
            return
    elif isinstance(diet,pd.DataFrame):
        print('Using custom diet!')
    else:
        logging.warn('Diet not used- please pass in valid DataFrame, local file, or resource file name!')
        return

    for m in model.reactions:
        if m.id.startswith('Diet_EX_'):
            m.bounds = (0,1000)

    d = get_adapted_diet(diet_df)
    logging.info('Adding %s diet to model...'%diet)
    added = []
    for m in model.reactions:
        f_id = m.id.split('[d]')[0].split('Diet_')[-1]
        if f_id in d.index:
            r = d.loc[f_id]
            m.bounds = (float(r.lb),float(r.ub))
            added.append(f_id)
    return added
