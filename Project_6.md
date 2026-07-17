# Cancer subtype classification from tumour gene expression

Cancer is not one disease but many, and two tumours that look identical under a microscope can behave very differently depending on which genes they express. Modern oncology increasingly uses these molecular 'subtypes' ( obtained from tumour gene-expression profiles) to guide diagnosis, prognosis, and treatment. Learning to classify tumours from expression data is a foundational skill in precision medicine. Enormous, openly accessible cancer-genomics resources make this possible without any special access or application. In this project you will take gene-expression data for a chosen cancer, use dimensionality reduction to reveal structure across samples, build a classifier for tumour subtype (or tumour versus normal tissue), and identify the genes that drive the prediction.

## Your task

1. Classify tumour subtype (or tumour vs normal, or a survival group) from expression profiles.
2. Identify the most informative genes and check biological plausibility.
3. Use PCA to reveal structure across samples.
4. Discuss what such a classifier would mean for diagnosis/stratification.

## Datasets

- **cBioPortal (TCGA and more)** — <https://www.cbioportal.org> — REST API. **Subset:** one cancer study.
- **NCBI GEO** — <https://www.ncbi.nlm.nih.gov/geo>. **Subset:** an expression series.
- **UCSC Xena** — <https://xena.ucsc.edu>. **Subset:** harmonised TCGA/GTEx; top-variance genes.

## Deliverables

Reproducible pipeline, PCA + subtype classifier with honest CV, an informative-gene list with biological interpretation, report with limitations (batch effects etc) and per-member roles.

## Stretch goals

Survival analysis; compare feature-selection strategies; train TCGA / test GEO cross-dataset validation.

---
*Work in a shared Git repo, split tasks via issues, and make your analysis fully
reproducible (pinned environment + scripted data download). Every member must
understand the whole project and state their role in the report.*
