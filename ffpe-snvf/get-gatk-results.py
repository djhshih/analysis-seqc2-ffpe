#!/usr/bin/env python
import polars as pl
from pathlib import Path
import vcformer as vcf
from tqdm import tqdm

repo_root = Path("..")

## Functions
def get_gatk_obmm(vcf_path: Path) -> pl.DataFrame:
	"""Reads VCF and returns the GATK Oritentation Bias mixture model's probability score and prediction"""
	gatk_obmm = (
		vcf.read_vcf_as_polars(vcf_path, info_fields=["ROQ"], sample_fields=[])
		.explode("filters", "alts")
		.rename({"filters" : "filter", "alts" : "alt"})
		.rename(lambda x: x.lower())
		.select("chrom", "pos", "ref", "alt", "filter", "roq")
		.with_columns(obmm_prob = 1 - (10 ** ( -pl.col("roq") / 10 )))
	)

	return gatk_obmm


def process_dataset(dataset: str, variant_set: str) -> None:
	"""Read the variant and filter columns from VCF and write 
	a table for evaluating the gatk orientation bias mixture model.

	Parameters:
		dataset: Dataset identifier (e.g. 'FFX').
		variant_set: Variant set folder name (e.g. 'mutect2-tn_filtered_pass-orientation').
	"""
	print(f"Processing Dataset : {dataset} | Variant Set : {variant_set}")

	outdir_root = (repo_root / "ffpe-snvf" / dataset / variant_set / "gatk-obmm")
	vcf_paths = list((repo_root / "vcf" / dataset / variant_set).glob("*/*.vcf.gz"))

	for path in tqdm(vcf_paths):
		sample_name = path.parent.name
		outdir = (outdir_root / sample_name)
		outdir.mkdir(exist_ok=True, parents=True)

		gatk_res = get_gatk_obmm(path)
		gatk_res.write_csv((outdir / sample_name).with_suffix(".gatk-obmm.tsv"), separator="\t")

## Get GATK OBMM results
process_dataset(dataset = "FFX", variant_set = "mutect2-tn_filtered_pass-orientation")
process_dataset(dataset = "FFG", variant_set = "mutect2-tn_filtered_pass-orientation")
