# Gut microbiome signatures of disease (metagenomics)

The trillions of microbes living in the human gut (the gut microbiome) shape health in ways we are only beginning to understand, with links to malnutrition, diarrhoeal disease, and metabolic and immune conditions of particular importance in Africa. Metagenomic sequencing lets us read which microbes are present in a sample and in what proportions, opening the door to finding microbial 'signatures' of disease. In this project you will take open microbiome data, ready-made taxonomic abundance tables, or raw sequencing reads for the more ambitious, and ask whether a person's disease state can be read from their microbial community. You will compare the diversity and composition of microbiomes between disease and control groups, identify the microbes that best distinguish them, and build a classifier, while grappling with the statistical quirks of microbiome data such as its compositional nature.

## Your task

1. Compare gut microbial composition between a disease group and controls.
2. Identify the most discriminative taxa and check biological plausibility.
3. Classify disease state from a taxonomic abundance profile.
4. Describe diversity patterns (alpha/beta) distinguishing the groups.

## Datasets

- **curatedMetagenomicData** — <https://waldronlab.io/curatedMetagenomicData> — Bioconductor package. **Subset:** one study's abundance table (recommended core).
- **NCBI SRA** — <https://www.ncbi.nlm.nih.gov/sra>. **Subset:** one small study of raw reads.
- **MGnify (EBI)** — <https://www.ebi.ac.uk/metagenomics>. **Subset:** pre-computed analyses.

## Deliverables

Reproducible workflow, alpha/beta-diversity analysis, a discriminative-taxa list, a disease-state classifier with honest evaluation, report with limitations (batch effects, compositionality) + per-member roles.

## Stretch goals

Process raw reads end-to-end on HPC; cross-study validation; functional (pathway) profiling.

---
*Work in a shared Git repo, split tasks via issues, and make your analysis fully
reproducible (pinned environment + scripted data download). Every member must
understand the whole project and state their role in the report.*
