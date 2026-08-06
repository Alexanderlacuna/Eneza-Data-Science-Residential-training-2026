# Eneza — Thyroid PTC subtype classification from RNA-seq

A reproducible machine-learning pipeline that classifies three histologic
subtypes of thyroid papillary carcinoma (PTC) from TCGA RNA-seq data:
classical/usual, follicular variant, and tall-cell variant.

This is a research prototype for molecular triage, not a clinical diagnostic
device.


## TLDR 

docker compose up --build

What happens:

1. Docker builds an image with Python 3.11 and the pinned dependencies.
2. The current project directory is mounted into `/app` inside the container.
3. The container runs `python scripts/run.py`.
4. All outputs are written back to the host in `./output/`.




## Quick start

### Docker (recommended for full reproducibility)

```bash
cd ~/Eneza-Data-Science-Residential-training-2026/Project6
docker compose up --build
```

Outputs are written to `output/`.

### Local Python

```bash
pip install -r requirements.txt
python scripts/run.py
```

Runtime: ~3–5 minutes on a laptop with 4 cores. See `REPRODUCIBILITY.md` for
Docker details and CI examples.

## Project structure

```text
.
├── data/                  # Raw TCGA data files
├── scripts/
│   ├── run.py             # Main pipeline
│   └── __init__.py
├── docs/                  # Reports and design notes
├── reports/               # PDF/HTML deliverables
├── output/                # Generated figures, tables, gene lists
├── requirements.txt       # Flexible dependencies
├── requirements.lock      # Pinned dependencies for Docker
├── Dockerfile
├── docker-compose.yml
├── Makefile
├── REPRODUCIBILITY.md
└── README.md
```

If you use version control, add `data/` and `output/` to `.gitignore`.

## Data

- `data/data_mrna_seq_v2_rsem.txt` — RSEM expression matrix (20,531 genes × 500 samples)
- `data/data_clinical_sample.txt` — sample labels (`TUMOR_TYPE`)
- `data/data_clinical_patient.txt` — patient survival data

After filtering, 491 samples remain: 355 classical, 102 follicular variant, 36
tall-cell variant.

## Pipeline

`scripts/run.py` runs an end-to-end workflow:

1. Clean and log2-transform expression data.
2. Select the top 5,000 variable genes.
3. Stratified 80/20 train/test split.
4. 5-fold stratified cross-validation with a scikit-learn pipeline:
   `SelectKBest(ANOVA) → StandardScaler → LogisticRegression(class_weight="balanced")`.
5. Hold-out evaluation, top-gene extraction, feature-selector and classifier
   benchmarks, and survival analysis.

Feature selection is fit inside each CV fold to avoid leakage. All random
operations use `random_state=42` and the full configuration is saved to
`output/tables/results_report.json`.

## Results

Default model (logistic regression, 100 ANOVA-selected genes):

| Metric | CV mean | Hold-out test |
|--------|---------:|---------------:|
| Balanced accuracy | 0.607 | 0.530 |
| Macro F1 | 0.550 | 0.503 |
| ROC-AUC OvR weighted | 0.747 | 0.752 |
| PR-AUC OvR weighted | — | 0.732 |

Performance is moderate, which is expected because the three PTC variants are
histologically and molecularly similar. The model is best used as a triage
prototype to flag samples for further review.

Full results, feature-selector/classifier comparisons, and survival analysis
are in `docs/PROJECT_REPORT.md` and `docs/anova_vs_deseq.org`.

## Outputs

```text
output/
├── figures/           # PCA, confusion matrix, top genes, KM curves, benchmarks
├── tables/            # JSON report, classification report, benchmark tables
└── gene_lists/        # Top predictive genes with biological notes
```

## Reproducibility

- Docker image uses Python 3.11 and `requirements.lock`.
- All random seeds are fixed.
- The `Config` dataclass is serialised into the JSON report.

To change the model or feature count, edit the `Config` object in
`scripts/run.py` or call `run_workflow()` with a custom `Config`.

## Limitations

- Single-cohort analysis (TCGA THCA only); external validation is needed.
- Small tall-cell class (36 samples) limits learning for that subtype.
- Global top-5,000 variance filter is unsupervised; move it inside CV for a
  production pipeline.
- Survival analysis is under-powered because PTC has very few death events.

## Citation and data provenance

- TCGA data: The Cancer Genome Atlas (TCGA) Thyroid Carcinoma (THCA) cohort,
  accessed via cBioPortal and TCGA clinical sample sheets.
- Design decisions: `docs/architecture_design.org`.

