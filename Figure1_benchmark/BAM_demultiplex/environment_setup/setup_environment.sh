#!/bin/bash
set -euo pipefail

echo "=== 1. Package Manager Setup ==="
# If conda/mamba is not available in the system, install micromamba in the ~/.local/bin directory of the current user.
PKG="conda"
RUNNER="conda run"

if command -v mamba >/dev/null 2>&1; then 
    PKG="mamba"
    RUNNER="mamba run"
elif ! command -v conda >/dev/null 2>&1; then
    echo "Conda/Mamba not found. Installing Micromamba locally..."
    BIN="$HOME/.local/bin"
    mkdir -p "$BIN"
    curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C "$BIN" --strip-components=1 bin/micromamba
    PKG="$BIN/micromamba"
    RUNNER="$BIN/micromamba run"
    export PATH="$BIN:$PATH"
fi

echo "=== 2. Creating Single Unified Conda Environment ==="
# Place all the dependencies in an environment called "hla_splitter_env"
ENV_NAME="hla_splitter_env"

"$PKG" create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    python=3.10 \
    pip \
    samtools \
    bcftools \
    htslib \
    zlib \
    bedtools \
    pigz \
    cellsnp-lite \
    vireosnp \
    sinto \
    subset-bam \
    kallisto=0.46.1 \
    'bustools>=0.43.1' \
    popscle \
    apptainer \
    arcas-hla \
    pandas

echo "=== Setup Complete! ==="
echo "To activate your environment, use:"
echo "conda activate $ENV_NAME   # (If using micromamba: micromamba activate $ENV_NAME)"
