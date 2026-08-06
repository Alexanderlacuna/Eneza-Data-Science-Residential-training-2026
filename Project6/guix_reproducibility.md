# Guix reproducibility alternative

This document explains how to run the Eneza pipeline reproducibly with **Guix**
instead of Docker. Guix is useful when Docker is unavailable, when you need
source-level provenance, or when you are working on an HPC cluster.

## Why Guix?

- **Functional package management:** Guix records the exact dependency graph and
  build recipes, not just installed binaries.
- **Declarative manifests:** A `manifest.scm` + `channels.scm` describes the
  whole environment and can be recreated with one command.
- **No root required:** `guix shell --container` can isolate the runtime
  without root privileges or a Docker daemon.
- **Rollback:** Guix profiles are immutable and can be rolled back instantly.
- **HPC friendly:** Guix is often available on academic clusters.

## Quick start

1. Make sure Guix is installed and the project data files are in `data/`.
2. Create a Guix manifest that includes Python and the required packages. A
   minimal starting point:

   ```scheme
   ;; manifest.scm
   (specifications->manifest
     (list "python@3.11"
           "python-numpy"
           "python-pandas"
           "python-scikit-learn"
           "python-matplotlib"
           "python-seaborn"
           "python-lifelines"))
   ```

3. Run the pipeline inside a Guix container:

   ```bash
   guix shell --container --network \
     --manifest=manifest.scm \
     --share=. \
     -- python scripts/run.py
   ```

   Outputs will be written to `output/`.

## Note on package versions

Guix package versions may differ from the versions pinned in `requirements.lock`
for Docker. For strict bit-for-bit reproducibility, use the Docker image. Use
Guix when you need source-level provenance, cluster compatibility, or cannot use
Docker.

## More information

- [Guix manual](https://guix.gnu.org/manual/)
- [Guix containers](https://guix.gnu.org/manual/en/html_node/Invoking-guix-shell.html)
