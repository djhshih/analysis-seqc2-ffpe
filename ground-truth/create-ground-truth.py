#!/usr/bin/env python
import polars as pl
import glob
import os
import logging
import sys

repo_root = ".."

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

	For each FFPE sample, compares variants to matched frozen samples from the same case,
	marks variants present in any frozen sample as 'truth', and writes the result to disk.

	Args:
		dataset: Dataset identifier.
		variant_set: Variant set identifier.
		ff_paths: Path to matched Fresh Frozen VCF.
		vcf_suffix: Suffix for VCF files (default: '.vcf.gz').
	"""
	logging.info(f"Processing dataset: {dataset}, variant set: {variant_set}")
	

	# Find all VCF file paths
	ffpe_paths = sorted(glob.glob(f"{repo_root}/vcf/{dataset}/{variant_set}/*/*{vcf_suffix}"))
	
	logging.info(f"Found {len(ffpe_paths)} FFPE VCF files")

	## Declare Fresh Frozen variant set name
	wes_ff_set = variant_set.replace("orientation", "orientation-exome") if "exome" not in variant_set else variant_set
	wgs_ff_set = variant_set.replace("-exome", "")

	## Glob Fresh Frozen VCF paths
	ff_paths = sorted(
		glob.glob(f"{repo_root}/vcf/WES/{wes_ff_set}/*/*.vcf.gz") +
		glob.glob(f"{repo_root}/vcf/WGS/{wgs_ff_set}/*/*.vcf.gz")
	)

	for i, path in enumerate(ffpe_paths, start=1):
		# Extract FFPE sample name and case ID
		ffpe_sample_name = path.split("/")[-2]

		logging.info(f"{i}. Processing FFPE sample: {ffpe_sample_name}")

		# Read and filter FFPE variants
		ffpe = (
			read_variants(path)
			.pipe(snv_filter)
			.select(["chrom", "pos", "ref", "alt"])
		)

		logging.info(f"\t{len(ff_paths)} matched FF samples found.")

		for j, ff_sample_path in enumerate(ff_paths, start=1):

			ff_sample = ff_sample_path.split("/")[-2]
			ff_col_name = f"in_ff_{j}"

			# Read frozen sample, filter variants
			logging.info(f"\t\tComparing with frozen sample {j}: {ff_sample}")

			ff = (
				read_variants(ff_sample_path)
				.pipe(snv_filter)
				.with_columns(pl.lit(True).alias(ff_col_name))
			)

			# Join FFPE and frozen sample on variant coordinates
			ffpe = (
				ffpe
				.join(ff, how="left", on=["chrom", "pos", "ref", "alt"])
				.with_columns(pl.col(ff_col_name).fill_null(False))
			)

		# Mark variants present in any frozen sample as "truth"
		ff_cols = [col for col in ffpe.columns if "in_ff_" in col]
		ffpe = ffpe.with_columns(pl.any_horizontal(ff_cols).alias("truth"))

		# Write ground-truth table to output directory
		outdir = f"{repo_root}/ground-truth/{dataset}/{variant_set}/{ffpe_sample_name}"
		os.makedirs(outdir, exist_ok=True)
		outpath = f"{outdir}/{ffpe_sample_name}.ground-truth.tsv"
		ffpe.write_csv(outpath, separator="\t")
		logging.info(f"\tGround-truth written to: {outpath}")
		

## Get Ground Truth
process_dataset(
    dataset="FFX",
    variant_set="mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist-clonal"
)

process_dataset(
    dataset="FFG",
    variant_set="mutect2-tn_filtered_pass-orientation-dp20-blacklist-clonal"
)


