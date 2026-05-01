# Secondary Structure Analysis Pipeline

This repository contains a Snakemake-based pipeline for analyzing protein secondary structure and computing backbone angles from PDB files.

## Overview

The pipeline performs the following steps:

* Unzips `.pdb.gz` files
* Runs STRIDE to extract secondary structure
* Processes STRIDE output with contextual information
* Computes backbone angles
* Generates plots for analysis

## Directory Structure

.
├── secondary_structure_analysis.smk
├── run.sh
├── scripts/
│   ├── stride_output_w_context.py
│   ├── compute_angles_from_tsv.py
│   └── plot_angles.py

## Requirements

### System Dependencies

* Python (>= 3.8)
* STRIDE
* gzip

### Python Dependencies

* pandas
* numpy
* matplotlib

### Snakemake

Install via conda:
conda install -c bioconda snakemake

## Input

* Input files should be in `.pdb.gz` format
* Update the PDB directory path in the Snakefile or config before running

## Running the Pipeline

### Run full pipeline

snakemake --snakefile secondary_structure_analysis.smk -j <num_cores>

### Run for a specific file

snakemake --snakefile secondary_structure_analysis.smk results/angles/<pdb_id>.txt -j 1

### Using run script

bash run.sh

## Output

* Final angle files: `results/angles/`
* Intermediate files: generated during execution
* Plots: saved as image files

## Notes

* Ensure STRIDE is available in PATH
* Use multiple cores for faster execution
* Large datasets may take significant time

## Author

Developed as part of coursework/research work in protein structure analysis.
