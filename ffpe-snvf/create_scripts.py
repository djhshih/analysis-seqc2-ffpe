import os
import glob
import polars as pl
import sys
import logging

## Local dependencies
repo_root = ".."
sys.path.append(f"{repo_root}/common-ffpe-snvf/python")
from common import return_path_if_exists

## Set up logging configuration
logging.basicConfig(
	stream=sys.stdout, 
	level=logging.INFO, 
	format='%(message)s', # %(asctime)s %(levelname)s:
	force=True
)

## Functions
def create_filtering_scripts(models: list, vcf_path: str, bam_path: str, filtered_outdir: str, sample_name: str, ref_path: str) -> None:
	for model in models:

		outdir = f"script_{model}"
		os.makedirs(outdir, exist_ok=True)
		template_dir = snvf_templates.get(model, None)
		filename = f"{outdir}/{model}_{sample_name}.sh"
	
		if model == "ideafix":
			content = [
				f"#!/usr/bin/env Rscript\n",
				f"Rscript {template_dir} --vcf '{vcf_path}' --ref '{ref_path}' --outdir '{filtered_outdir}/{model}/{sample_name}'\n"
			]
		
		elif model == "ffpolish":
			content = [
				"#!/usr/bin/env bash \n",
				f"ffpolish filter -o {filtered_outdir}/{model}/{sample_name} -p {sample_name} {ref_path} {vcf_path} {bam_path} \n"
			]

		else:
			content = [
				"#!/bin/bash\n",
				f"bash {template_dir} '{bam_path}' '{vcf_path}' '{filtered_outdir}/{model}/{sample_name}' '{ref_path}'\n"
			]
		
		with open(filename, "w") as f:
			f.writelines(content)
			
		logging.info(f"\t\t- Created {model.upper()} filtering script")

def process_dataset(dataset: str, variant_set: str) -> None:
	vcf_paths = sorted([os.path.abspath(path) for path in glob.glob(f"{repo_root}/vcf/{dataset}/{variant_set}/*/*.vcf.gz")])
	filtered_outdir = os.path.abspath(f"{repo_root}/ffpe-snvf/{dataset}/{variant_set}")
	bam_dir = f"{repo_root}/data/{dataset}/bam"

	logging.info(f"Found {len(vcf_paths)} VCFs to process...")

	for vcf_path in vcf_paths:

		logging.info(f"\tProcessing: {vcf_path}")
		sample_name = os.path.basename(vcf_path).split(".")[0]
		
		bam_path = return_path_if_exists(f"{bam_dir}/{sample_name}.bwa.dedup.bam", abs=True)
		create_filtering_scripts(models, vcf_path, bam_path, filtered_outdir, sample_name, ref_path)

#--------------------

## Setup
ref_path = return_path_if_exists(f"{repo_root}/data/seqc2-reference-genome/GRCh38/GRCh38.d1.vd1.fa", abs=True)

models = ["mobsnvf", "vafsnvf", "sobdetector", "ideafix", "ffpolish"]

snvf_templates = {
    "ideafix" : return_path_if_exists(f"{repo_root}/common-ffpe-snvf/R/ideafix.R", abs=True),
    "mobsnvf" : return_path_if_exists(f"{repo_root}/common-ffpe-snvf/templates/ffpe-snvf/mobsnvf.sh.template", abs=True),
    "sobdetector" : return_path_if_exists(f"{repo_root}/common-ffpe-snvf/templates/ffpe-snvf/sobdetector.sh.template", abs=True),
    "vafsnvf" : return_path_if_exists(f"{repo_root}/common-ffpe-snvf/templates/ffpe-snvf/vafsnvf.sh.template", abs=True)
}

## Create filtering batch scripts
process_dataset("FFX", "mutect2-tn_filtered_pass-orientation")
process_dataset("FFG", "mutect2-tn_filtered_pass-orientation")

