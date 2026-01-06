#!/usr/bin/env python
import os
import glob
import polars as pl

# %%
def return_path_if_exists(path: str, abs: bool =False) -> str:
	if abs:
		path = os.path.abspath(path)
	if os.path.exists(path):
		return path
	else:
		raise FileNotFoundError(f"Path {path} does not exist.")

# %%
ref_path = return_path_if_exists("../data/seqc2-reference-genome/GRCh38/GRCh38.d1.vd1.fa", abs=True)
vcf_paths = [os.path.abspath(path) for path in sorted(glob.glob("../vcf/mutect2-matched-normal_pass-orientation-dp-filtered/FF*/*/*.vcf"))]

filtered_outdir = return_path_if_exists(".", abs=True)

model = "ffpolish"


# %%

for vcf_path in vcf_paths:
	
	sample_name = vcf_path.split("/")[-2]
	dataset = vcf_path.split("/")[-3]
	vcf_type = vcf_path.split("/")[-4]
	
	bam_path = return_path_if_exists(glob.glob(f"../data/bam/{dataset}/{sample_name}*.bam")[0], abs=True)

	content = f"#!/usr/bin/env bash\nffpolish filter -o {filtered_outdir}/{vcf_type}/{dataset}/ffpolish/{sample_name} -p {sample_name} {ref_path} {vcf_path} {bam_path}\n"
	
	batch_scripts_outdir = "script_ffpolish"
	os.makedirs(batch_scripts_outdir, exist_ok=True)
	filename = f"{batch_scripts_outdir}/{model}_{sample_name}.sh"
	
	with open(filename, "w") as f:
		f.write(content)

		print(f"Created {model} filtering script for {vcf_path}")

