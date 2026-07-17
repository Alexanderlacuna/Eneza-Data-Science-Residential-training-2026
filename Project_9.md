# Resolving cryptic coral diversity in Kenyan reefs using low-coverage whole-genome sequencing

## Abstract

Coral species identification underpins reef monitoring and conservation, yet coral morphology is highly plastic, and many reef-building genera harbour cryptic, genetically distinct lineages that are grouped under a single name. The Western Indian Ocean (WIO) remains poorly characterised genomically. This project provides a compact real-world dataset: 60 coral colonies from four reef sites near Shimoni, Kenya, spanning five genera (Acropora, Pocillopora, Stylophora, Porites, Millepora), sequenced using low-coverage whole-genome sequencing (lcWGS). From the same reads, three independent data layers can be extracted: mitochondrial barcodes (COX1, mtORF e.t.c), ultraconserved elements (UCE), and genome-wide nuclear SNPs, enabling a cross-validated approach to species identity. The overarching objective is species delimitation: determining whether morphologically recognised taxa correspond to evolutionarily distinct lineages, and identifying cryptic diversity that standard surveys would miss.

## Research Objectives

1.	Investigate anomalous samples phylogenetic placement and gene-flow tests.
2.	Compile a provenance-tracked reference panel recording type status (holotype / topotype / vetted non-type) for every reference, and re-run placements against type-anchored sequences to harden species names.
3.	Profile Symbiodiniaceae communities (ITS2) per host lineage to test whether cryptic host lineages carry distinct symbionts.
4.	Apply formal species-delimitation methods (e.g. mPTP, BPP, or SNAPP) to quantify whether the divergent lineages represent a distinct species.

## Expected output
Participants will construct a reproducible  workflow pipeline that takes raw reads through quality control, mapping, organellar assembly, mitochondrial gene trees, UCE phylogenomics, and low-coverage population structure (ANGSD, NGSadmix), and tests for gene flow. 

    version-controlled Snakemake/nextfow pipeline, annotated trees and population-structure plots, a type-status reference table, symbiont profiles per lineage, delimitation results, and a short group report proposing provisional (cf./aff.) species names.

## Dataset

## Useful resources