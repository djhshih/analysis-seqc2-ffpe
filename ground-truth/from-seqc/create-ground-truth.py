#!/usr/bin/env python
import polars as pl
import glob
import os
import logging
import sys

repo_root = "../.."

## Local Dependencies
sys.path.append(f"{repo_root}/common-ffpe-snvf/python")
from common import read_variants, snv_filter


## Set up logging configuration
logging.basicConfig(
	stream=sys.stdout, 
	level=logging.INFO, 
	format='%(message)s', # %(asctime)s %(levelname)s:
	force=True
)

## Functions
def process_dataset(dataset: str, variant_set: str, vcf_suffix: str = ".vcf.gz") -> None:
	"""
	Process a dataset to generate ground-truth variant tables for FFPE samples.

	For each FFPE sample, compares variants to matched high confidence variantss from the same case,
	marks variants present in any high confidence variants as 'truth', and writes the result to disk.

	Args:
		dataset: Dataset identifier.
		variant_set: Variant set identifier.
		ff_paths: Path to matched Fresh Frozen VCF.
		vcf_suffix: Suffix for VCF files (default: '.vcf.gz').
	"""
	logging.info(f"Processing dataset: {dataset}, variant set: {variant_set}")
	

	# Find all VCF file paths
	ffpe_paths = sorted(glob.glob(f"{repo_root}/vcf/{dataset}/{variant_set}/*/*{vcf_suffix}"))
	
	logging.info(f"\tFound {len(ffpe_paths)} FFPE VCF files")

	## ground truth path
	high_conf_path = f"{repo_root}/ground-truth/from-seqc/seqc2_ground-truth_snv_pass-high-conf.tsv"

	for i, path in enumerate(ffpe_paths, start=1):
		# Extract FFPE sample name and case ID
		ffpe_sample_name = path.split("/")[-2]

		logging.info(f"\t{i}. Processing FFPE sample: {ffpe_sample_name}")

		# Read and filter FFPE variants
		ffpe = (
			read_variants(path)
			.pipe(snv_filter)
			.select(["chrom", "pos", "ref", "alt"])
		)

		# Read high confidence variants, and mark them true
		logging.info(f"\t\tComparing with high confidence variant list: {high_conf_path}")
		high_conf = pl.read_csv(high_conf_path, separator="\t").drop("filter").with_columns(pl.lit(True).alias("truth"))

		# Join FFPE and high confidence variants on variant coordinates
		ffpe = (
			ffpe
			.join(high_conf, how="left", on=["chrom", "pos", "ref", "alt"])
			.with_columns(pl.col("truth").fill_null(False))
		)

		# Write ground-truth table to output directory
		outdir = f"{repo_root}/ground-truth/from-seqc/{dataset}/{variant_set}/{ffpe_sample_name}"
		os.makedirs(outdir, exist_ok=True)
		outpath = f"{outdir}/{ffpe_sample_name}.ground-truth.tsv"
		ffpe.write_csv(outpath, separator="\t")
		logging.info(f"\t\tGround-truth written to: {outpath}\n")
		

## Get Ground Truth

### FFX
process_dataset(
    dataset="FFX",
    variant_set="mutect2-tn_filtered_pass-orientation"
)

process_dataset(
    dataset="FFX",
    variant_set="mutect2-tn_filtered_pass-orientation-dp20"
)

process_dataset(
    dataset="FFX",
    variant_set="mutect2-tn_filtered_pass-orientation-dp20-blacklist"
)

process_dataset(
    dataset="FFX",
    variant_set="mutect2-tn_filtered_pass-orientation-exome"
)

process_dataset(
    dataset="FFX",
    variant_set="mutect2-tn_filtered_pass-orientation-exome-dp20"
)

process_dataset(
    dataset="FFX",
    variant_set="mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist"
)

process_dataset(
    dataset="FFX",
    variant_set="mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist-clonal"
)

### FFG
process_dataset(
    dataset="FFG",
    variant_set="mutect2-tn_filtered_pass-orientation"
)

process_dataset(
    dataset="FFG",
    variant_set="mutect2-tn_filtered_pass-orientation-dp20"
)

process_dataset(
    dataset="FFG",
    variant_set="mutect2-tn_filtered_pass-orientation-dp20-blacklist"
)

process_dataset(
    dataset="FFG",
    variant_set="mutect2-tn_filtered_pass-orientation-dp20-blacklist-clonal"
)
