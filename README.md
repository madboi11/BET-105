# Secondary Structure Analysis Pipeline

This repository contains a Snakemake-based pipeline for analyzing protein secondary structure and computing backbone angles from PDB files. The workflow processes compressed PDB files, extracts structural features, computes geometric angles, and generates distribution plots.

## Overview

The pipeline performs the following steps:

* Unzips `.pdb.gz` files
* Runs STRIDE to assign secondary structure
* Formats STRIDE output with context
* Computes backbone angles
* Generates KDE-based plots

## Directory Structure

```
.
├── secondary_structure_analysis.smk
├── run.sh
├── scripts/
│   ├── stride_output_w_context.py
│   ├── compute_angles_from_tsv.py
│   └── plot_angles.py
│
├── results/                 # created automatically
│   ├── angles/              # final angle files (.txt)
│   ├── plots/               # generated plots (.png)
│   └── valid_pdbs.txt       # list of valid PDBs used
│
├── unzipped_pdbs            # intermediate unzipped pdbs (temporary files)
├── stride_out/              # intermediate STRIDE outputs
├── stride_w_context/        # intermediate context TSV files

```

## Input

* Input files must be in `.pdb.gz` format
* Update the `pdb_dir` path when running Snakemake

---

## Running the Pipeline

Once inside the correct environment and current working directory,

```bash
bash run.sh path/to/data
```

## Output

* `results/angles/` → computed angle values
* `results/plots/` → final KDE plots
* `results/valid_pdbs.txt` → list of PDBs that contributed to plots
