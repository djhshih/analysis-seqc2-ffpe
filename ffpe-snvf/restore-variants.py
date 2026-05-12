#!/usr/bin/env python
import polars as pl
import glob
import os
import sys
from tqdm import tqdm

repo_root = ".."
## Local Dependencies
sys.path.append(f"{repo_root}/common-ffpe-snvf/python")
from common import read_variants, snv_filter, ct_filter

# This script aims to restore the variants omitted by ffpolish and ideafix during evaluation 
# and assigns appropiate score based the conditions under which they were omitted

## Functions
def restore_ffpolish(df: pl.DataFrame) -> pl.DataFrame:
	return (
		df.with_columns(
			pl.when(pl.col("filter") != "PASS")
			.then(pl.lit(0).alias("score"))
			.otherwise(pl.col("score"))
		)
		.with_columns(
			pl.when(pl.col("filter") != "PASS")
			.then(pl.lit(0).alias("ffpolish_pred"))
			.otherwise(pl.col("ffpolish_pred"))
		)
	)

def restore_ideafix(df: pl.DataFrame) -> pl.DataFrame:
	return (
		df
		.with_columns(
			pl.when(pl.col("filter") != "PASS")
			.then(pl.lit(1).alias("deam_score"))
			.otherwise(pl.col("deam_score")),

			pl.when(pl.col("filter") != "PASS")
			.then(pl.lit("non-PASS").alias("deamination"))
			.otherwise(pl.col("deamination")),
		)
		.with_columns(
			pl.when(pl.col("deam_score").is_null())
			.then(pl.lit(0).alias("deam_score"))
			.otherwise(pl.col("deam_score")),

			pl.when(pl.col("deamination").is_null())
			.then(pl.lit("VAF>0.3").alias("deamination"))
			.otherwise(pl.col("deamination")),
		)
	)

def process_dataset(dataset, variant_set):

	print(f'Processing Dataset: {dataset} | Variant Set: {variant_set}')

	to_fix_paths = sorted(
		glob.glob(f"{repo_root}/ffpe-snvf/{dataset}/{variant_set}/ffpolish/*/*.ffpolish.tsv") +
		glob.glob(f"{repo_root}/ffpe-snvf/{dataset}/{variant_set}/ideafix/*/*.ideafix-xgboost.tsv")
	)

	for path in tqdm(to_fix_paths):
		sample_name = path.split("/")[-2]
		model = path.split("/")[-3]

		snvf = pl.read_csv(path, separator="\t")

		vcf_path = f"{repo_root}/vcf//{dataset}/{variant_set}/{sample_name}/{sample_name}.vcf.gz"
		vcf = read_variants(vcf_path, columns=["#CHROM", "POS", "REF", "ALT", "FILTER"]).pipe(ct_filter)

		if model == "ffpolish":
			restored_snvf = vcf.pipe(snv_filter).join(snvf, on=["chrom", "pos", "ref", "alt"], how="left").pipe(restore_ffpolish)

		if "ideafix" in model:
			restored_snvf  = vcf.pipe(ct_filter).join(snvf, on=["chrom", "pos", "ref", "alt"], how="left").pipe(restore_ideafix)

		unrestored_path = ".".join(path.split(".")[:-1] + ["incomplete"] + [path.split(".")[-1]])
		os.rename(path, unrestored_path)

		restored_snvf.write_csv(path, separator="\t")


## Restore the variants
process_dataset("FFX", "mutect2-tn_filtered_pass-orientation")

process_dataset("FFG", "mutect2-tn_filtered_pass-orientation")


