# Single-cell transcriptomics of disease

Every tissue in the body is a mixture of many different cell types, each switching different genes on and off. Traditional 'bulk' sequencing measures the average signal across all of them at once, hiding the very differences between cells that matter most in disease.Single-cell RNA sequencing (scRNA-seq) measures gene expression in thousands of individual cells simultaneously, letting researchers discover cell types, track how they shift in disease, and pinpoint the marker genes that define each cell state. Think of it like bulk RNASeq looks at a fruit cocktail smooth as is while ScRNA-seq looks at the individual fruits in the smoothy. Large, openly available scRNA-seq datasets now exist for many diseases and tissues. In this project you will take one such open, disease-relevant dataset and run a complete exploratory pipeline (quality control to remove low-quality cells, clustering cells into groups, identifying the marker genes that label each group, and compare healthy and diseased states. The end product is a map of which cells are present in a tissue and how they change in disease.

## Your task

1. Justify QC thresholds and report cells/genes retained.
2. Cluster cells and identify marker genes for each cluster.
3. Assign biological cell-type identities.
4. Find differentially expressed genes between diseased and healthy state.

## Datasets

- **CELLxGENE Discover** — <https://cellxgene.cziscience.com> via the `cellxgene-census` client or direct H5AD (no request). **Subset:** one small disease-relevant dataset (a few thousand cells).
- **NCBI GEO** — <https://www.ncbi.nlm.nih.gov/geo>. **Subset:** a scRNA-seq series.
- **Human Cell Atlas** — <https://data.humancellatlas.org>.

## Deliverables

Reproducible notebook (raw matrix → QC → clustering → annotation), UMAP plots, cluster marker tables, an annotated cell-type figure, a report with a short biological interpretation and per-member roles.

## Stretch goals

Differential expression between conditions; compare two clustering resolutions; automated vs manual annotation.

---
*Work in a shared Git repo, split tasks via issues, and make your analysis fully
reproducible (pinned environment + scripted data download). Every member must
understand the whole project and state their role in the report.*
