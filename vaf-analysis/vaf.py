#!/usr/bin/env python
import vcformer as vcf
import polars as pl
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import glob

# %% [markdown]
# ## Functions

def transpose(df: pl.DataFrame) -> pl.DataFrame:
	temp = df.transpose()
	stats = df.transpose(column_names = temp.row(0))[1:]
	return stats

def process_vcf(vcf_path: str, vcf_type: str, sample_name: str, normal_name: str = "", snv_only = True) -> pl.DataFrame:
	# only get the Allelic Fractions field i.e the proportion of read supporting the ALT allele
	df = vcf.read_vcf_as_polars(vcf_path, sample_fields=['AF', "AD"], info_fields=[])
	df = (
		df
		.rename({"alts" : "alt"})
		.rename(lambda x : x.replace(f"{sample_name}", "tumor").replace(f"{normal_name}", "normal"))
		.with_columns(pl.col(["alt", "filters", "tumor.AF"]).explode())
		# .with_columns(pl.lit(sample_name).alias("sample_name"))
	)

	if snv_only:
		df = df.filter((pl.col("ref").str.len_chars() == 1)  & (pl.col("alt").str.len_chars() == 1))

	# Calculate VAF from AD as AF is based on a probabilistic model 
	df = df.with_columns((pl.col("tumor.AD").list.get(1) / (pl.col("tumor.AD").list.get(0) + pl.col("tumor.AD").list.get(1))).alias("tumor.AD_derived_VAF"))
 		
	df = df.with_columns(pl.lit(sample_name).alias("sample_name"))
 
	if vcf_type == "matched-normal":
		cols = ['sample_name', 'chrom', 'pos', 'ref', 'alt', 'tumor.AF', 'tumor.AD', 'normal.AF', 'normal.AD', 'tumor.AD_derived_VAF']
	else:
		cols = ['sample_name', 'chrom', 'pos', 'ref', 'alt', 'tumor.AF',  'tumor.AD', 'tumor.AD_derived_VAF']
  
	df = df.select(cols)

	return df

def perform_analysis(df: pl.DataFrame, vcf_type: str, sample_name: str, outdir: str, write: bool = True) -> tuple[pl.DataFrame, pl.DataFrame]:

	if write:
		df.write_parquet(f"{outdir}/{sample_name}.{vcf_type}.variants-af-table.parquet")

	stats_af = df[f"tumor.AF"].describe()
	if write:
		stats_af.write_csv(f"{outdir}/{sample_name}.{vcf_type}.AF.stats.tsv", separator="\t")
	stats_af = transpose(stats_af).with_columns(pl.lit(sample_name).alias("sample_name"))
	stats_af = stats_af.select(['sample_name', 'count', 'null_count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'])

	stats_vaf = df[f"tumor.AD_derived_VAF"].describe()
	if write:
		stats_vaf.write_csv(f"{outdir}/{sample_name}.{vcf_type}.AD_derived_VAF.stats.tsv", separator="\t")
	stats_vaf = transpose(stats_vaf).with_columns(pl.lit(sample_name).alias("sample_name"))
	stats_vaf = stats_vaf.select(['sample_name', 'count', 'null_count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'])
	
	if write:
		sns.displot(data=df, x=f"tumor.AF", bins = 100).savefig(f"{outdir}/{sample_name}.{vcf_type}.AF.pdf")
		plt.close()
		sns.displot(data=df, x=f"tumor.AD_derived_VAF", bins = 100).savefig(f"{outdir}/{sample_name}.{vcf_type}.AD_derived_VAF.pdf")
		plt.close()
	else:
		sns.displot(data=df, x=f"tumor.AF", bins = 100)
		sns.displot(data=df, x=f"tumor.AD_derived_VAF", bins = 100)

	return stats_af, stats_vaf

# %% [markdown]
# ### VAF Analysis
# 
# - AF is the allele fractions field in the mutect2 VCF which applies a probabilistic model to obtain VAF
# - AD is the allelic depth field which shows reads supporting each allele. VAF is also derived from here: `alt AD / (ref AD + alt AD)`
# 
# Individual result sets using each analyses are produced

dataset_normals = {
	"FFG" : "WGS_IL_N_1",
	"WGS" : "WGS_IL_N_1",
	"FFX" : "WES_IL_N_2",
	"WES" : "WES_IL_N_2"
}

datasets = ["WES", "FFX", "WGS", "FFG"]

for dataset in datasets:
	
	print(f"Processing dataset: {dataset}")

	normal_name = dataset_normals[dataset]
	outdir_root = f"./{dataset}"

	matched_vcf_paths = glob.glob(f"../vcf/mutect2-matched-normal_pass-orientation-filtered/{dataset}/*/*.vcf")
	tumor_only_vcf_paths = glob.glob(f"../vcf/mutect2-tumor-only_pass-orientation-filtered/{dataset}/*/*.vcf")

	tumor_only_vcf_table = pl.DataFrame({
		"sample_name" : [os.path.basename(path).split(".")[0] for path in tumor_only_vcf_paths],
		"tumor_only_vcf_path" : tumor_only_vcf_paths
	})

	matched_normal_vcf_table = pl.DataFrame({
		"sample_name" : [os.path.basename(path).split(".")[0] for path in matched_vcf_paths],
		"matched_normal_vcf_path" : matched_vcf_paths
	})

	vcf_table = matched_normal_vcf_table.join(tumor_only_vcf_table, on="sample_name").sort("sample_name")
	vcf_table

	matched_normal_combined = []
	tumor_only_combined = []
	matched_stats_af = []
	matched_stats_vaf = []
	tumor_only_stats_af = []
	tumor_only_stats_vaf = []

	for i, sample_name in enumerate(vcf_table["sample_name"]):

		outdir = f"{outdir_root}/{sample_name}"
		os.makedirs(outdir, exist_ok=True)

		print(f"{i+1}. Processing {sample_name}")
	
		matched = process_vcf(vcf_table[i, "matched_normal_vcf_path"], "matched-normal", sample_name, normal_name)
		stats_af, stats_vaf = perform_analysis(matched, "matched-normal", sample_name, outdir)
		matched_stats_af.append(stats_af)
		matched_stats_vaf.append(stats_vaf)
		
		tumor_only = process_vcf(vcf_table[i, "tumor_only_vcf_path"], "tumor-only", sample_name, normal_name)
		stats_af, stats_vaf = perform_analysis(tumor_only, "tumor-only", sample_name, outdir)
		tumor_only_stats_af.append(stats_af)
		tumor_only_stats_vaf.append(stats_vaf)

		matched_normal_combined.append(matched)
		tumor_only_combined.append(tumor_only)

	matched_normal_combined = pl.concat(matched_normal_combined, how="diagonal_relaxed")
	tumor_only_combined = pl.concat(tumor_only_combined, how="diagonal_relaxed")

	matched_stats_af = pl.concat(matched_stats_af)
	tumor_only_stats_af = pl.concat(tumor_only_stats_af)

	matched_stats_vaf = pl.concat(matched_stats_vaf)
	tumor_only_stats_vaf = pl.concat(tumor_only_stats_vaf)

	print("Performing combined analysis...")

	matched_stats_af.write_csv(f"{outdir_root}/per-{dataset.lower()}-samples.matched-normal.AF.stats.tsv", separator="\t")
	tumor_only_stats_af.write_csv(f"{outdir_root}/per-{dataset.lower()}-samples.tumor-only.AF.stats.tsv", separator="\t")

	matched_stats_vaf.write_csv(f"{outdir_root}/per-{dataset.lower()}-samples.matched-normal.AD_derived_VAF.stats.tsv", separator="\t")
	tumor_only_stats_vaf.write_csv(f"{outdir_root}/per-{dataset.lower()}-samples.tumor-only.AD_derived_VAF.stats.tsv", separator="\t")

	_, _ = perform_analysis(matched_normal_combined, "matched-normal", f"all-{dataset.lower()}-samples", outdir_root)
	_, _ = perform_analysis(tumor_only_combined, "tumor-only", f"all-{dataset.lower()}-samples", outdir_root)

