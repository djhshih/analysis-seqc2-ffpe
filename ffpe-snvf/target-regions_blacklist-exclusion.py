#!/usr/bin/env python
import polars as pl
import os
import glob
from tqdm import tqdm 

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

def get_filtering_summary(
    sample_name: str, 
    model: str, 
    orig_var_set: str, 
    new_var_set: str, 
    snvf: pl.DataFrame, 
    target_vars: pl.DataFrame, 
    filtered_snvf: pl.DataFrame
) -> dict:

    n_orig = snvf.height
    n_target = target_vars.height
    n_filtered = filtered_snvf.height

    summary = {
        "original_var_set" : orig_var_set,
        "new_var_set" : new_var_set,
        "sample_name" : sample_name,
        "model" : model,
        "n_var_original" : n_orig,
        "n_var_target_vcf" : n_target,
        "n_var_filtered_snvf": n_filtered,
        "n_var_removed" : n_orig - n_filtered,
        "pct_removed" : ((n_orig - n_filtered) / n_orig) * 100
    }

    return summary

# %% [markdown]
# ## Setup

repo_root = ".."

wgs_vcf_dir = f"{repo_root}/vcf/mutect2-matched-normal_hc-target_pass-orientation-dp-blacklist-filtered"
wes_vcf_dir = f"{repo_root}/vcf/mutect2-matched-normal_exome-hc-target_pass-orientation-dp-blacklist-filtered"

# %% [markdown]
# ## WES SNVF Filtering

## Subset_WES Data
wes_ffpe_snvf = (
	glob.glob("mutect2-matched-normal_pass-orientation-dp-filtered_micr1234/FFX/*/*/*.tsv") +
	glob.glob("mutect2-matched-normal_pass-orientation-dp-filtered_micr1234/FFX/*/*/*.snv")
)

new_var_set = "mutect2-matched-normal_exome-hc-target_pass-orientation-dp-blacklist-filtered"
wes_filtering_summary = []

for path in tqdm(wes_ffpe_snvf):
	dataset = path.split("/")[-4]
	model = path.split("/")[-3]
	sample_name = path.split("/")[-2]
	variant_set = path.split("/")[-5]
	fname = os.path.basename(path)

	snvf = pl.read_csv(path, separator="\t", infer_schema_length=1000)
	
	target_vars_path = f"{wes_vcf_dir}/{dataset}/{sample_name}/{sample_name}.bwa.dedup.vcf"
	target_vars = read_variants(target_vars_path)

	filtered_snvf = snvf.join(target_vars, on = ["chrom", "pos", "ref", "alt"], how="semi")

	filtered_snvf_outdir = f"{new_var_set}/{dataset}/{model}/{sample_name}"
	os.makedirs(filtered_snvf_outdir, exist_ok=True)

	sample_filtering_summary = get_filtering_summary(
		sample_name,
		model,
		variant_set,
		new_var_set,
		snvf,
		target_vars,
		filtered_snvf
	)

	wes_filtering_summary.append(sample_filtering_summary)

	filtered_snvf.write_csv(f"{filtered_snvf_outdir}/{fname}", separator="\t")

pl.DataFrame(wes_filtering_summary).write_csv(f"{new_var_set}/target-regions_blacklist-exclusion_filtering-summary.tsv", separator="\t")



# %% [markdown]
# ## WGS SNVF Filtering

## Subset_wgs Data
wgs_ffpe_snvf = (
	glob.glob("mutect2-matched-normal_pass-orientation-dp-filtered/FFG/*/*/*.ffpolish.tsv") +
	glob.glob("mutect2-matched-normal_pass-orientation-dp-filtered/FFG/*/*/*.gatk-obmm.tsv") +
	glob.glob("mutect2-matched-normal_pass-orientation-dp-filtered/FFG/*/*/*.ideafix-xgboost.tsv") +
	glob.glob("mutect2-matched-normal_pass-orientation-dp-filtered/FFG/*/*/*.ideafix-rf.tsv") +
	glob.glob("mutect2-matched-normal_pass-orientation-dp-filtered/FFG/*/*/*.mobsnvf.snv") +
	glob.glob("mutect2-matched-normal_pass-orientation-dp-filtered/FFG/*/*/*.sobdetector.snv") +
	glob.glob("mutect2-matched-normal_pass-orientation-dp-filtered/FFG/*/*/*.vafsnvf.snv")
)

new_var_set = "mutect2-matched-normal_hc-target_pass-orientation-dp-blacklist-filtered"
wgs_filtering_summary = []

for path in tqdm(wgs_ffpe_snvf):
	dataset = path.split("/")[-4]
	model = path.split("/")[-3]
	sample_name = path.split("/")[-2]
	variant_set = path.split("/")[-5]
	fname = os.path.basename(path)

	snvf = pl.read_csv(path, separator="\t", infer_schema_length=1000)
	
	target_vars_path = f"{wgs_vcf_dir}/{dataset}/{sample_name}/{sample_name}.bwa.dedup.vcf"
	target_vars = read_variants(target_vars_path)

	filtered_snvf = snvf.join(target_vars, on = ["chrom", "pos", "ref", "alt"], how="semi")

	filtered_snvf_outdir = f"{new_var_set}/{dataset}/{model}/{sample_name}"
	os.makedirs(filtered_snvf_outdir, exist_ok=True)

	sample_filtering_summary = get_filtering_summary(
		sample_name,
		model,
		variant_set,
		new_var_set,
		snvf,
		target_vars,
		filtered_snvf
	)

	wgs_filtering_summary.append(sample_filtering_summary)

	filtered_snvf.write_csv(f"{filtered_snvf_outdir}/{fname}", separator="\t")

pl.DataFrame(wgs_filtering_summary).write_csv(f"{new_var_set}/target-regions_blacklist-exclusion_filtering-summary.tsv", separator="\t")




