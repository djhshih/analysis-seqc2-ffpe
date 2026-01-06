#!/usr/bin/env python
import polars as pl
import os
import glob
from tqdm import tqdm

def read_variants(path:str, columns: list = ["#CHROM", "POS", "REF", "ALT", "FILTER"]) -> pl.DataFrame:
	variants = (
		pl.read_csv(path, separator="\t", comment_prefix="##", infer_schema_length=1000, columns=columns)
		.rename(lambda x: x.lstrip("#").lower())
		.with_columns(pl.col("alt").str.split(","))
		.explode("alt")
	)
	return variants

datasets = ["FFX", "FFG"]
variant_types = [
	"mutect2-matched-normal_pass-orientation-filtered",
	"mutect2-matched-normal_pass-orientation-dp-filtered", 
	"mutect2-matched-normal_pass-orientation-dp20-filtered",
	"mutect2-tumor-only_pass-orientation-filtered",
	"mutect2-tumor-only_pass-orientation-dp-filtered",
	"mutect2-tumor-only_pass-orientation-dp20-filtered"

]

for dataset in datasets:
	print(f"{dataset}:")
	for variant_type in variant_types:
		print(F"\t{variant_type}")
		vcf_paths = sorted(glob.glob(f"../vcf/{variant_type}/{dataset}/*/*.vcf"))

		for path in tqdm(vcf_paths):
			sample_name = path.split("/")[-2]

			variants = read_variants(path)

			outdir = f"{variant_type}/{dataset}/gatk-obmm/{sample_name}"
			os.makedirs(outdir, exist_ok=True)

			variants.write_csv(f"{outdir}/{sample_name}.tsv", separator="\t")



