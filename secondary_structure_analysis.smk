import glob
import os

PDB_DIR = config["pdb_dir"]
PDBS = [os.path.basename(f).replace(".pdb.gz", "") for f in glob.glob(os.path.join(PDB_DIR, "*.pdb.gz"))]

rule all:
    input:
        expand("results/angles/{pdb}.txt", pdb=PDBS)

rule unzip_pdb:
    input:
        pdb = os.path.join(PDB_DIR, "{pdb}.pdb.gz")
    output:
        unzip = temp("unzipped_pdbs/{pdb}.pdb")
    shell:
        "zcat {input.pdb} > {output.unzip}"

rule run_stride:
    input:
        unzip = "unzipped_pdbs/{pdb}.pdb"
    output:
        stride = "stride_out/{pdb}.stride"
    shell:
        "stride {input.unzip} > {output.stride} 2>/dev/null || touch {output.stride}"

rule format_stride:
    input:
        stride = "stride_out/{pdb}.stride"
    output:
        tsv = "stride_w_context/{pdb}_ARG.tsv"
    shell:
        "python scripts/stride_output_w_context.py {input.stride} {output.tsv} ARG"

rule compute_angles:
    input:
        tsv = "stride_w_context/{pdb}_ARG.tsv",
        pdb = "unzipped_pdbs/{pdb}.pdb"
    output:
        angles = "results/angles/{pdb}.txt"
    shell:
        "python scripts/compute_angles_from_tsv.py {input.tsv} {input.pdb} {output.angles}"
