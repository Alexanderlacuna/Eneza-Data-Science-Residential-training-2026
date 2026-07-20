# Predicting tuberculosis drug resistance from genomes

Tuberculosis is among the world's leading infectious killers and a major burden in Kenya. Diagnosing drug-resistant TB with traditional culture-based testing can take weeks. Whole-genome sequencing offers a much faster route: because particular mutations confer resistance to particular drugs, reading a bacterium's genome can predict its resistance profile in days rather than weeks. Large collections of *Mycobacterium tuberculosis* genomes, paired with resistance phenotypes and a WHO-curated catalogue of resistance mutations, are openly available. In this project you will profile the resistance mutations circulating in African isolates and build a model that predicts resistance from genomic features, benchmarking your predictions against the WHO catalogue.

## Your task

1. Catalogue which resistance mutations (*rpoB*, *katG*, *inhA*, *gyrA*, *pncA*) occur across African isolates and how common they are.
2. Predict resistance to a drug (e.g. rifampicin, isoniazid) from genomic features.
3. Compare your predictions to the WHO mutation catalogue.
4. Describe how resistance varies by country/lineage.

## Datasets

- **NCBI Pathogen Detection — *M. tuberculosis*** — <https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Mycobacterium_tuberculosis/>. **Subset:** a few hundred isolates from one region.
- **BV-BRC / PATRIC** — <https://www.bv-brc.org> — REST API. **Subset:** curated AMR phenotypes.
- **WHO TB mutation catalogue (2023)** — <https://www.who.int/publications>. **Subset:** reference truth set.
- **CRyPTIC data (EBI/figshare)** — <https://www.ebi.ac.uk>. **Subset:** genotype–phenotype.

## Deliverables

Reproducible pipeline (pinned env), resistance-mutation summary, a classifier evaluated against the WHO catalogue, report with limitations and per-member roles.

## Stretch goals

Multi-drug prediction; lineage calling; rules-based (catalogue) vs ML.

---
*Work in a shared Git repo, split tasks via issues, and make your analysis fully
reproducible (pinned environment + scripted data download). Every member must
understand the whole project and state their role in the report.*
