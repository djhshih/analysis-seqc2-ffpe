#!/usr/bin/env python
import os
import polars as pl

# %% [markdown]
# ## Functions

def read_variants(path:str, columns: list = ["#CHROM", "POS", "REF", "ALT"]) -> pl.DataFrame:
	"""
	Reads Variants from a VCF file into a Polars DataFrame.
	By default, only reads essential columns.
	Splits multi-allelic variants into separate rows.
	Can be extended to read additional columns as needed or
	read custom variant files with similar structure.
	
	:param path: Path to VCF file
	:type path: str
	:param columns: List of column names to read from the VCF file. Defaults to ["#CHROM", "POS", "REF", "ALT"].
	:type columns: list
	:return: DataFrame containing the variants
	:rtype: DataFrame
	"""
	variants = (
		pl.read_csv(path, separator="\t", comment_prefix="##", infer_schema_length=1000, columns=columns)
		.rename(lambda x: x.lstrip("#").lower())
		.with_columns(pl.col("alt").str.split(","))
		.explode("alt")
	)
	return variants

# %% [markdown]
# ## Main

repo_root = ".."
snv_gt_path = f"{repo_root}/data/release/latest/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
gt_outdir = f"{repo_root}/eval-seqc2-gt/ground-truth"
os.makedirs(gt_outdir, exist_ok=True)

snv_gt = read_variants(snv_gt_path, columns = ["#CHROM", "POS", "REF", "ALT", "FILTER"])
snv_gt_pass_hc = snv_gt.filter(pl.col("filter") == "PASS;HighConf")


snv_gt.write_csv(f"{gt_outdir}/seqc2_ground-truth_snv_pass-high-med-conf.tsv", separator="\t")
snv_gt_pass_hc.write_csv(f"{gt_outdir}/seqc2_ground-truth_snv_pass-high-conf.tsv", separator="\t")


