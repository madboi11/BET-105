#!/bin/bash
set -uo pipefail

if [ -z "${1:-}" ]; then
    echo "Usage: $0 /path/to/pdb_dir [cores]"
    exit 1
fi

PDB_DIR=$(realpath "$1")
CORES=${2:-$(nproc)}

mkdir -p unzipped_pdbs stride_out stride_w_context results/angles results/plots

snakemake \
    --snakefile secondary_structure_analysis.smk \
    --cores "$CORES" \
    --keep-going \
    --rerun-incomplete \
    --config pdb_dir="$PDB_DIR"

echo "Running plotting..."

python scripts/plot_angles.py results/angles results/plots/ss_profile_HHH_for_arg_with_valid_runs.png
