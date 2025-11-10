#!/usr/bin/env python
import os
import glob
import polars as pl


def assign_path_if_exist(path: str) -> str:
    if os.path.exists(path):
        return path
    else:
        raise FileNotFoundError(f"{path} does not exist")


ffpe_datasets = ["FFG", "FFX"]
models = ["mobsnvf", "vafsnvf", "sobdetector"]

## VCFs called using mutect2 with 1 matched normal
filtered_outdir_root = os.path.abspath("..")

outdir = "script"
os.makedirs(outdir, exist_ok=True)

vcf_paths = glob.glob("../../../vcf/mutect2-tumor-only_pass-orientation-filtered/*/*/*.vcf")

vcf_table = pl.DataFrame(
	{
		"sample_name" : [os.path.basename(path).split(".")[0] for path in vcf_paths],
		"vcf_path" : [os.path.abspath(path) for path in vcf_paths]
	}
)

bam_paths = glob.glob("../../../data/bam/*/*.bam")

bam_table = pl.DataFrame(
	{
		"sample_name" : [os.path.basename(path).split(".")[0] for path in bam_paths],
		"bam_path" : [os.path.abspath(path) for path in bam_paths]
	}
)

bam_vcf_table = bam_table.join(vcf_table, on="sample_name", how="inner")

# Removing non-ffpe datasets
bam_vcf_table = bam_vcf_table.filter(pl.col("sample_name").str.contains_any(ffpe_datasets))

for i in range(bam_vcf_table.shape[0]):
	
	bam_path = bam_vcf_table[i, "bam_path"]
	vcf_path = bam_vcf_table[i, "vcf_path"]
	sample_name = os.path.basename(vcf_path).replace(".vcf", "").replace(".gz", "")
	dset = bam_path.split("/")[-2]
	
	for model in models:
    
		filtered_outdir = f"{filtered_outdir_root}/{dset}"

		
		filter_script_path = assign_path_if_exist(f"../../{model}.sh.template")
		
		content = f"#!/bin/bash\nbash ../{filter_script_path} '{bam_path}' '{vcf_path}' '{filtered_outdir}'\n"
		filename = f"{outdir}/{model}_{sample_name}.sh"
		
		with open(filename, "w") as f:
			f.write(content)
			
		print(f"Created {model} filtering script for {vcf_path}")



