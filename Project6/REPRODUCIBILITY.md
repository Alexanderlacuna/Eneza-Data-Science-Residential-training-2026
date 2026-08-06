# Eneza group 6 — Reproducibility & Docker Deployment Guide

This guide explains how to run the Eneza Project 6  pipeline inside Docker so that the
analysis is fully reproducible across machines, operating systems, and time.

If you prefer using guix we have a section for you [./guix_reproducilibity.md]

---



## Why we chose Docker?

Reproducibility has three layers:

1. **Deterministic code** — seeded random numbers, immutable `Config`, and a
   single entry point (`scripts/run.py`).
2. **Pinned software environment** — exact Python interpreter and package
   versions.
3. **Fixed input data** — the same TCGA files mounted into the container.

Docker addresses layer 2 by bundling the Python version and all dependencies
into an image.  Docker Compose addresses layers 1 and 3 by mounting the project
folder and running the same command every time.

---

## Files added for reproducibility

| File | Purpose |
|------|---------|
| `Dockerfile` | Builds the reproducible Python runtime image |
| `docker-compose.yml` | One-command build + run service |
| `requirements.lock` | Exact versions of the direct Python dependencies |
| `.dockerignore` | Keeps large data files and outputs out of the image |
| `Makefile` | Short commands: `make build`, `make run`, `make shell`, `make freeze` |
| `docs/architecture_design.org` | Why we made  architectural decision  |

---

## Prerequisites

- [Docker Engine](https://docs.docker.com/engine/install/) 24.0+ (or Docker Desktop)
- [Docker Compose](https://docs.docker.com/compose/install/) v2+
- The Eneza repository cloned locally, **including the raw TCGA data files**
  you can fetch the file from 
  https://github.com/cBioPortal/datahub/tree/master/public/thca_tcga_pan_can_atlas_2018
  
  - `data/data_mrna_seq_v2_rsem.txt`
  - `data/data_clinical_sample.txt`
  - `data/data_clinical_patient.txt`

---

## Quick start

The fastest way to reproduce the full pipeline is with Docker Compose:

```bash
# From the repository root navigate to  i.e ./project6
docker compose up --build
```

What happens:

1. Docker builds an image with Python 3.11 and the pinned dependencies.
2. The current project directory is mounted into `/app` inside the container.
3. The container runs `python scripts/run.py`.
4. All outputs are written back to the host in `./output/`.

--- 

## Manual Docker commands

If you prefer plain Docker instead of Compose:

```bash
# 1. Build the image
docker build -t eneza .

# 2. Run the pipeline (mounts current directory into /app)
docker run --rm -v "$(pwd)":/app eneza

# 3. Inspect the outputs
ls output/tables
ls output/figures
ls output/gene_lists
```

---

## Interactive exploration

To drop into a shell inside the same reproducible environment:

```bash
# With Docker Compose
docker compose run --rm eneza bash

# With plain Docker
docker run --rm -it -v "$(pwd)":/app --entrypoint /bin/bash eneza
```

Inside the container you can run Python, inspect data, or re-run the pipeline
with a modified configuration.

---

## What gets produced

After a successful run the following artefacts are written to `output/` on the
host (the container writes through the volume mount):

```text
output/
├── figures/
│   ├── pca_subtype.png
│   ├── confusion_matrix.png
│   ├── top_genes.png
│   ├── feature_selector_comparison.png
│   ├── classifier_comparison.png
│   ├── km_annotated_subtype.png
│   └── km_predicted_subtype.png
├── tables/
│   ├── results_report.json
│   ├── classification_report.csv
│   ├── feature_selector_comparison.csv
│   ├── classifier_comparison.csv
│   ├── survival_logrank.csv
│   └── deseq_vs_anova.csv
└── gene_lists/
    └── top_predictive_genes.csv
```

`output/tables/results_report.json` contains the full `Config` object and all
scores, so every output directory is self-documenting.

---

## Pinning dependencies even harder

`requirements.lock` pins the **direct** dependencies that we install
explicitly (numpy, pandas, scikit-learn, matplotlib, seaborn, lifelines).
When pip installs them it also resolves transitive dependencies
(joblib, scipy, contourpy, etc.).

If you need a fully transitive lock file, run:

```bash
make freeze
```

This creates `requirements-freeze.txt` with every single package version that
is actually installed in the container.  You can then use it in the Dockerfile
instead of `requirements.lock` for maximum reproducibility.

---

## Changing the configuration

`scripts/run.py` defines a `Config` dataclass at the top of the script.  To run with
a different model or feature count, create a small wrapper script and mount it
into the container:

```python
# custom_run.py
from scripts.run import Config, run_workflow

cfg = Config(
    classifier="random_forest",
    selector="kbest",
    selected_features=250,
    cv_splits=5,
    random_state=42,
    output_dir="./output_random_forest",
)
run_workflow(cfg)
```

Run it with:

```bash
docker run --rm -v "$(pwd)":/app eneza python custom_run.py
```

Because `random_state=42` is fixed, the same config will produce identical
results every time.


## Troubleshooting

### `FileNotFoundError` for data files

The TCGA data files must be present in the `data/` directory before the
container starts.  They are not included in the Docker image.  Check that
`data/data_mrna_seq_v2_rsem.txt`, `data/data_clinical_sample.txt`, and
`data/data_clinical_patient.txt` exist.

you can fetch them from the TCGA repository:
https://github.com/cBioPortal/datahub/tree/master/public/thca_tcga_pan_can_atlas_2018

### Permission errors on output files

By default the container runs as root, so files in `output/` are owned by root.
On Linux you can change ownership afterwards:

```bash
sudo chown -R "$USER:$USER" output/
```

