# Eneza reproducible runtime
# --------------------------
# This image contains only the code and its pinned Python dependencies.
# TCGA data files are NOT baked into the image; they are mounted from the
# host at runtime.  This keeps the image small and guarantees that the same
# input files always produce the same outputs.

FROM python:3.11-slim

# Avoid interactive prompts and keep Python output unbuffered.
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

# Install the dependency lock first so Docker can cache this layer when only
# the code changes.
COPY requirements.lock /app/requirements.lock
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /app/requirements.lock

# Copy the Python package.  Data files are excluded by .dockerignore.
COPY scripts/ /app/scripts/

# Default behaviour: run the full reproducible pipeline.
ENTRYPOINT ["python", "scripts/run.py"]
CMD []
