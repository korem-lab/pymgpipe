import sys
import os
import gc
import tqdm
import pandas as pd
import logging
import matplotlib.pyplot as plt 
import numpy as np
import time
import skbio 
import seaborn as sns 
from cobra import Configuration
from scipy.spatial.distance import squareform, pdist
from pathlib import Path
from multiprocessing import Pool
from functools import partial
from .build import _build
from .diet import add_diet_to_model
from .io import load_cobra_model, write_lp_problem, write_cobra_model, suppress_stdout
from .utils import load_dataframe, remove_reverse_vars
from .coupling import add_coupling_constraints
from .metrics import compute_diversity_metrics

cobra_config = Configuration()
cobra_config.lower_bound = -1000
cobra_config.upper_bound = 1000


def build_models(
    coverage_file,
    taxa_dir,
    solver="gurobi",
    threads=int(os.cpu_count() / 2),
    samples=None,
    parallel=True,
    lp_type=".mps",
    cobra_type=".xml",
    out_dir="./",
    coupling_constraints=True,
    diet_fecal_compartments=False,
    remove_reverse_vars_from_lp=False,
    hard_remove=False,
    abundance_threshold=1e-6,
    diet=None,
    vaginal=False,
    essential_metabolites=None, 
    micronutrients=None,
    force_uptake=True,
    diet_threshold=0.8,
    compress=True,
    compute_metrics=True,
    force=False,
    sample_prefix='mc'
):
    """Build community COBRA models using mgpipe-like compartments and constraints.

    This function is pymgpipe's main model building function, and can be used to build models for either one or multiple samples

    Args:
        coverage_file (pandas.DataFrame | str): Abundance matrix with taxa as rows and samples as columns
        taxa_dir (str): Directory containing individual strain/species taxa models (file names corresponding to index of coverage matrix)
        out_dir (str): Directory to save output of this function (models, LP problems, etc.), defaults to cwd
        diet_fecal_compartments (bool): Build models with mgpipe's diet/fecal compartmentalization, defaults to False
        coupling_constraints (bool): Add mgpipe's abundance coupling constraints to internal reactions, defaults to True
        remove_reverse_vars_from_lp (bool): Remove reverse variables that are added into COBRApy models
        hard_remove (bool): If set to True, variables will be removed entirely from models, which can be a time consuming process. Having it set to False will `remove` them by setting their bounds to 0.
        diet (str): Name of diet or file you want to impose on models (can be one diet or personalized)
        vaginal (bool): Impose vaginal diet constraints to support vaginal diet
        essential_metabolites (list): List of metabolites to include as essential within the diet
        micronutrients (list): List of metabolites to include as micronutrients within diet
        force_uptake (bool): Force minimum uptake of nutrients from diet
        diet_threshold (float): Threshold on imposed diet uptake, i.e. lb < metabolite < (threshold * lb)
        solver (str): LP solver (gurobi or cplex) used to solve models, defaults to gurobi
        parallel (bool): Samples will be built in parallel if set to True
        threads (int): Number of threads to use if building in parallel
        lp_type (str): File type for LP problem (either .mps or .lp), defaults to .mps
        cobra_type (str): File type for COBRA model (.xml, .mat, .json), defaults to .xml
        compress (bool): Models and LP problems will be saved as compressed files if set to True, defaults to True
        compute_metrics (bool): Compute diversity metrics for built models, defaults to True

    Notes:
        COBRA models written to *out_dir/models/*\n
        LP problems written to *out_dir/problems/*

    """
     
    start = time.time()
    cobra_config.solver = solver

    gc.disable()
    Path(out_dir).mkdir(exist_ok=True)

    taxa_dir = taxa_dir + "/" if taxa_dir[-1] != "/" else taxa_dir
    out_dir = out_dir + "/" if out_dir[-1] != "/" else out_dir

    model_dir = out_dir + "models/"
    problem_dir = out_dir + "problems/"
    Path(model_dir).mkdir(exist_ok=True)
    Path(problem_dir).mkdir(exist_ok=True)

    formatted = _format_coverage_file(coverage_file, out_dir, sample_prefix)
    samples_to_run = list(formatted.columns)
    taxa = list(formatted.index)

    print(
        "Found coverage file with %s samples and %s unique taxa"
        % (len(samples_to_run), len(taxa))
    )

    if samples is not None:
        samples_to_run = samples if isinstance(samples, list) else [samples]

    if len(samples_to_run) <= 1:
        parallel = False

    threads = os.cpu_count() - 1 if threads == -1 else threads
    threads = min(threads, len(samples_to_run))

    print("Building %s models..." % len(samples_to_run))
    print("\n----------------Parameters----------------")
    print("Diet/fecal compartments- %s" % str(diet_fecal_compartments).upper())
    print("Coupling constraints- %s" % str(coupling_constraints).upper())
    print("Parallel- %s" % str(parallel).upper())
    print("Threads- %s" % str(threads).upper())
    print("Solver- %s" % solver.upper())
    print("LP type- %s" % str(lp_type.split(".")[1]).upper())
    print("COBRA type- %s" % str(cobra_type.split(".")[1]).upper())
    print("compress- %s" % str(compress).upper())
    print("Output directory- %s" % str(out_dir).upper())

    _func = partial(
        _build_single_model,
        formatted,
        taxa_dir,
        solver,
        model_dir,
        problem_dir,
        lp_type,
        cobra_type,
        coupling_constraints,
        diet_fecal_compartments,
        remove_reverse_vars_from_lp,
        hard_remove,
        diet,
        vaginal,
        essential_metabolites,
        micronutrients,
        force_uptake,
        diet_threshold,
        abundance_threshold,
        compress,
        compute_metrics,
        force
    )

    if parallel:
        p = Pool(threads, initializer=_mute)
        p.daemon = False

        metrics = list(
            tqdm.tqdm(p.imap(_func, samples_to_run), total=len(samples_to_run))
        )
        p.close()
        p.join()
    else:
        metrics = tqdm.tqdm(list(map(_func, samples_to_run)), total=len(samples_to_run))

    if compute_metrics:
        try:
            abundance_df = pd.DataFrame({m['sample']:m['reaction_abundance'] for m in metrics})
            abundance_df.to_csv(out_dir+'reaction_abundance.csv')

            abundance_df[abundance_df > 0] = 1
            abundance_df.to_csv(out_dir+'reaction_content.csv')

            unique_reactions = [(len(m['taxa']),len(m['unique_reactions'])) for m in metrics]

            plt.scatter(*zip(*unique_reactions))
            plt.xlabel('# Taxa', fontsize=12)
            plt.ylabel('# Unique Reactions', fontsize=12)
            plt.xticks(np.arange(0,max([r[0] for r in unique_reactions])+1,1))
            plt.title('Metabolic Diversity')
            plt.savefig(out_dir+'metabolic_diversity.png')

            plt.clf()

            # PCoA plot of reaction presence 
            abundance_df = abundance_df.fillna(0).T
            distance_matrix = squareform(pdist(abundance_df, 'braycurtis'))
            my_pcoa = skbio.stats.ordination.pcoa(distance_matrix)

            pcoa_points = my_pcoa.samples
            pcoa_points.set_index(abundance_df.index,inplace=True)

            ax = sns.scatterplot(data=pcoa_points, x='PC1', y='PC2', hue=pcoa_points.index, ax=plt.gca())
            str1 = 'PC1 ' + str(round(my_pcoa.proportion_explained[0], 3) * 100) + '% of explained variance'
            str2 = 'PC2 ' + str(round(my_pcoa.proportion_explained[1], 3) * 100) + '% of explained variance'
            ax.set(xlabel=str1, ylabel=str2, title='PCoA of reaction presence (braycurtis)')
            plt.legend([],[], frameon=False)

            for line in range(0,pcoa_points.shape[0]):
                ax.text(pcoa_points['PC1'][line]+0.001, pcoa_points['PC2'][line]+0.001, 
                pcoa_points.index[line], horizontalalignment='left', 
                size='small', color='black', weight='light')

            plt.savefig(out_dir+'reaction_pcoa.png')


        except Exception as e:
            logging.warn('Ran into problem while saving metrics...skipping this step!\n%s'%e)
            pass
        
    print("-------------------------------------------------------")
    print("Finished building %s models and associated LP problems!" % len(metrics))
    print('Process took %s minutes to run...'%round((time.time()-start)/60,3))


def _build_single_model(
    coverage_df,
    taxa_dir,
    solver,
    model_dir,
    problem_dir,
    lp_type,
    cobra_type,
    coupling_constraints,
    diet_fecal_compartments,
    remove_reverse_vars_from_lp,
    hard_remove,
    diet,
    vaginal,
    essential_metabolites, 
    micronutrients,
    force_uptake,
    diet_threshold,
    abundance_threshold,
    compress,
    compute_metrics,
    force,
    sample_label,
):
    model_out = (
        model_dir + "%s.%s" % (sample_label, cobra_type.split(".")[1])
        if not compress
        else model_dir + sample_label + ".xml.gz"
    )
    lp_out = problem_dir + "%s.%s" % (sample_label, lp_type.split(".")[1])

    pymgpipe_model = None
    metrics = None

    if not force and os.path.exists(model_out):
        pymgpipe_model = load_cobra_model(model_out)
    else:
        with suppress_stdout():
            pymgpipe_model = _build(
                sample=sample_label,
                abundances=coverage_df,
                taxa_directory=taxa_dir,
                threshold=abundance_threshold,
                diet_fecal_compartments=diet_fecal_compartments,
                solver=solver
            )
        force = True 
        write_cobra_model(pymgpipe_model, model_out)

    if compute_metrics:
        metrics = compute_diversity_metrics(pymgpipe_model)
        if metrics is None or len(metrics)==0:
            logging.warning('Unable to compute diversity metrics for %s'%pymgpipe_model.name)

    if force or (not os.path.exists(lp_out) and not os.path.exists(lp_out+'.gz') and not os.path.exists(lp_out+'.7z')):
        # ----- START OPTLANG MODIFICATIONS -----
        if remove_reverse_vars_from_lp:
            try:
                logging.info('Removing variables!')
                remove_reverse_vars(pymgpipe_model,hard_remove)
            except Exception:
                logging.warning('Failed to remove reverse variables!')

        if diet is not None:
            add_diet_to_model(pymgpipe_model, diet, force_uptake, essential_metabolites, micronutrients, vaginal, diet_threshold)

        if coupling_constraints:
            try:
                add_coupling_constraints(pymgpipe_model)
            except Exception:
                logging.warning("Failed to add coupling constraints!")

        write_lp_problem(pymgpipe_model, out_file=lp_out, compress=compress, force=True)
    else:
        logging.info('Skipping %s because LP problem already exists!'%sample_label)
    del pymgpipe_model
    gc.collect()
    return metrics


def _format_coverage_file(coverage_file, out_dir = './', sample_prefix = None):
    coverage = load_dataframe(coverage_file)

    if sample_prefix is None:
        return coverage
    
    conversion_file_path = out_dir + "sample_label_conversion.csv"
    try:
        sample_conversion_dict = pd.read_csv(conversion_file_path, index_col=0).iloc[:, 0].to_dict()
        if set(sample_conversion_dict.keys()) != set(coverage.columns):
            raise Exception(
                "Provided label conversion file %s does not provide labels for all samples!"
                % conversion_file_path
        )
    except:
        sample_conversion_dict = {
            v: sample_prefix + str(i + 1) for i, v in enumerate(sorted(coverage.columns))
        }
        pd.DataFrame({"conversion": sample_conversion_dict}).to_csv(conversion_file_path)
    
    return coverage.rename(columns=sample_conversion_dict)

    # melted = pd.melt(
    #     coverage,
    #     value_vars=coverage.columns,
    #     ignore_index=False,
    #     var_name="sample_id",
    #     value_name="abundance",
    # )
    # melted["strain"] = melted.index
    # melted["file"] = melted.strain.apply(
    #     lambda x: existing_taxa_files[x] if x in existing_taxa_files else None
    # )
    # melted["original_id"] = melted.sample_id.apply(lambda x: conversion_t[x])

    # missing_taxa = melted[melted.file.isna()].index.unique()
    # if len(missing_taxa) > 0:
    #     melted.drop(missing_taxa, axis="index", inplace=True)
    #     print("Removed %s missing taxa from coverage file-" % len(missing_taxa))
    #     print(missing_taxa)

    # melted = melted.loc[melted.abundance != 0]
    # melted.abundance = melted.abundance.div(
    #     melted.groupby(["sample_id"])["abundance"].transform("sum")
    # )

    # melted.reset_index(inplace=True)
    # melted.rename({"index": "id"}, axis="columns", inplace=True)

    # return melted

def _mute():
    sys.stdout = open(os.devnull, "w")