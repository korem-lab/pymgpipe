---
title: 'pymgpipe: microbiome metabolic modeling in Python'
tags:
  - Python
  - metabolic modeling
  - flux balance analysis
  - microbiome
  - cobra
  - mgpipe
authors:
  - name: Yoli Meydan
    orcid: 0009-0003-4597-3340
    equal-contrib: true
    affiliation: 1
  - name: Federico Bladini
    orcid: 0000-0001-9079-8869
    equal-contrib: true 
    affiliation: 1
  - name: Tal Korem
    orcid: 0000-0002-0609-0858
    corresponding: true 
    affiliation: "1, 2"
affiliations:
 - name: Program for Mathematical Genomics, Department of Systems Biology, Columbia University Irving Medical Center, New York, NY, USA
   index: 1
 - name: Department of Obstetrics and Gynecology, Columbia University Irving Medical Center, New York, NY, USA
   index: 2
date: 16 May 2023
bibliography: paper.bib
---

# Summary

Microbially-produced metabolites and microbiome metabolism in general are strongly linked to ecosystem-level phenotypes, including the health of the human host [1], [2]. To aid in the study of microbial metabolism from observational, human-derived data, a variety of computational methods that predict microbial community metabolic output from taxonomic abundances have been developed [3]–[6]. Several of these methods rely on community-scale metabolic models, which are mechanistic, knowledge-based models that enable the formulation and in silico testing of biological hypotheses regarding the metabolism of microbial communities [4], [5]. Community-scale models primarily use Flux Balance Analysis, a modeling technique that infers the metabolic fluxes in a system by optimizing an objective function, typically growth rate, subject to an assumption of a steady state and constraints imposed by the metabolic reactions present in the system [7]. These metabolic reactions are obtained from genome-scale metabolic networks (GEMs), knowledge-based computational models encompassing the known biochemical   reactions present within an organism [8]. In recent years, curated GEMs for thousands of human-associated microbial organisms have become increasingly available, allowing for a more in-depth exploration of the human microbiome [9]–[11].  In recent years, several community-scale metabolic modeling methods specifically tailored to the human microbiome have emerged, such as MICOM and mgPipe [4], [5].