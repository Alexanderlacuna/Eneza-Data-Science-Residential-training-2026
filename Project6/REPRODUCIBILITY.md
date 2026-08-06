# Eneza — Reproducibility & Deployment Guide

This guide describes the two reproducible execution contexts supported for the
Eneza thyroid-cancer subtype pipeline.

Reproducibility has three layers:

1. **Deterministic code** — seeded random numbers, immutable `Config`, and a
   single entry point (`scripts/run.py`).
2. **Pinned software environment** — exact Python interpreter and package
   versions.
3. **Fixed input data** — the same TCGA files mounted into the runtime.

## Supported contexts

| Context | Best for | How it pins the environment |
|---------|----------|----------------------------|
| **Docker** (default) | Laptops, shared servers, CI/CD, reviewers who just want a single command. | `Dockerfile` + `requirements.lock` bundle Python 3.11 and all dependencies into an image. |
| **Guix** (alternative) | Academic clusters, HPC, environments where Docker is unavailable or where source-level provenance is required. | `channels.scm` + `manifest.scm` describe the full package dependency graph and build recipes. |

---

## Docker (default)

Docker is the primary reproducibility context. It is the easiest way to share
an identical environment across different machines.

### Files added for Docker

| File | Purpose |
|------|---------|
| `Dockerfile` | Builds the reproducible Python runtime image |
| `docker-compose.yml` | One-command build + run service |
| `requirements.lock` | Exact versions of the direct Python dependencies |
| `.dockerignore` | Keeps large data files and outputs out of the image |
| `Makefile` | Short commands: `make build`, `make run`, `make shell`, `make freeze` |

### Prerequisites

- [Docker Engine](https://docs.docker.com/engine/install/) 24.0+ (or Docker Desktop)
- [Docker Compose](https://docs.docker.com/compose/install/) v2+
- The raw TCGA data files in `data/`:
  - `data/data_mrna_seq_v2_rsem.txt`
  - `data/data_clinical_sample.txt`
  - `data/data_clinical_patient.txt`

You can fetch the TCGA THCA Pan-Cancer Atlas files from the
[cBioPortal datahub](https://github.com/cBioPortal/datahub/tree/master/public/thca_tcga_pan_can_atlas_2018).

### Quick start with Docker Compose

```bash
# From the repository root
docker compose up --build
```

What happens:

1. Docker builds an image with Python 3.11 and the pinned dependencies.
2. The current project directory is mounted into `/app` inside the container.
3. The container runs `python scripts/run.py`.
4. All outputs are written back to the host in `./output/`.

Expected runtime: **3–5 minutes** on a modern laptop.

### Manual Docker commands

```bash
# Build the image
docker build -t eneza .

# Run the pipeline
docker run --rm -v "$(pwd)":/app eneza

# Inspect outputs
ls output/tables
ls output/figures
ls output/gene_lists
```

### Interactive shell

```bash
docker compose run --rm eneza bash
```

or

```bash
docker run --rm -it -v "$(pwd)":/app --entrypoint /bin/bash eneza
```

### Pinning dependencies even harder

`requirements.lock` pins the direct dependencies. To create a fully transitive
lock file, run:

```bash
make freeze
```

This produces `requirements-freeze.txt`, which you can use in place of
`requirements.lock` for maximum reproducibility.

### Changing the configuration

`scripts/run.py` defines an immutable `Config` dataclass at the top of the file.
To run a custom configuration, create a wrapper script and mount it into the
container:

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

```bash
docker run --rm -v "$(pwd)":/app eneza python custom_run.py
```

---

## Guix (alternative)

For environments where Docker is unavailable or undesirable, the project can
also be run reproducibly with **Guix**, a functional package manager. Guix gives
source-level provenance, declarative manifests, and containerised execution
without root privileges.

See the full guide: **[guix_reproducibility.md](guix_reproducibility.md)**.

In short: **Docker is the easiest portable default; Guix is the strongest
provenance-and-auditability alternative**, especially for HPC and shared
clusters.

---

## What gets produced

After a successful run the following artefacts are written to `output/`:

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

## Troubleshooting

### `FileNotFoundError` for data files

The TCGA data files must be present in `data/` before the container or Guix
environment starts. Check that the following files exist:

- `data/data_mrna_seq_v2_rsem.txt`
- `data/data_clinical_sample.txt`
- `data/data_clinical_patient.txt`

You can download them from the [cBioPortal datahub](https://github.com/cBioPortal/datahub/tree/master/public/thca_tcga_pan_can_atlas_2018).

### Permission errors on output files

By default the Docker container runs as root, so `output/` files may be owned
by root. On Linux you can fix this with:

```bash
sudo chown -R "$USER:$USER" output/
```

Guix containers do not have this issue because they run under your own UID.

---

## Summary

- Run with Docker: `docker compose up --build`.
- Run with Guix: see [`guix_reproducibility.md`](guix_reproducibility.md).
- Either context produces identical outputs from the same `data/` files and the
  same `scripts/run.py` entry point.
