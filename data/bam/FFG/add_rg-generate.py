#!/usr/bin/env python
import os
import glob
import polars as pl

## Additional info can be obtained through the annotations. However, this is not necessary

# annot = pl.read_csv("../../../annot/sample_annotation.tsv", separator="\t").filter(pl.col("study_alias") == "FFG")
# annot

# ffg_meta = (
# 		pl.read_excel("../../../annot/seqc2_dataset-annotations.xlsx", sheet_name="FFG")
# 		.rename(lambda x : x.lower().replace(" ", "_"))
# 	)
# ffg_meta

bam_paths = glob.glob("*.bam")
bam_paths

script_header="""#!/usr/bin/env bash

# This script is used to add missing @RG to the BAM files

set -euo pipefail

mkdir -p no_rg_bam
mkdir -p rg_fix_temp
"""

script_outdir = "batch_scripts"
os.makedirs(script_outdir, exist_ok=True)

for path in bam_paths:
    file_name = os.path.basename(path)
    sample_name = file_name.split(".")[0]
    
    ## Library (LB), Platform Unit (PU), Center Name (CN) maybe be specified but it is not necessary for out use case
    ## These can be obtained from the annotation table
    rg=f"ID:{sample_name}\\tLB:na\\tPL:ILLUMINA\\tSM:{sample_name}\\tPU:na\\tCN:na\\tDS:{file_name}"
    
    rg_replace_line = f"samtools addreplacerg -@ 4 -r $'{rg}' -o rg_fix_temp/{file_name} {path}"
    index_line = f"samtools index -@ 4 -b -o rg_fix_temp/{file_name.replace(".bam", ".bai")} rg_fix_temp/{file_name}"
    move_old = f"mv {path} no_rg_bam/"
    move_new = f"mv rg_fix_temp/{file_name} rg_fix_temp/{file_name.replace(".bam", ".bai")} ./"
    
    full_script = f"{script_header}\n{rg_replace_line}\n{index_line}\n{move_old}\n{move_new}\n"
    
    script_path = f"{script_outdir}/{sample_name}-add_rg.sh"
    
    with open(script_path, "w") as file:
        file.write(full_script)



