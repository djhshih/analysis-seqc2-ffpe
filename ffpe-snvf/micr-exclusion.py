#!/usr/bin/env python
import os
import glob
import polars as pl
from tqdm import tqdm

## Functions
def check_res_existence(path: str, abs = False, halt: bool = False) -> str | None:
	if os.path.exists(path):
		return os.path.abspath(path) if abs else path
	else:
		if halt:
			raise FileNotFoundError(f"{path} does not exist.")
		else:
			print(f"'{path}' does not exist. skipping...")
			return None


def filter_micr_artifacts(ffpe_filter: str, dataset: str, new_set_name: str, source_set_name: str, sample: str, micr_artifacts: pl.DataFrame, res_suffix: str = ".snv") -> pl.DataFrame:

	print(f"\tFiltering MICR artifacts from {ffpe_filter.upper()} results...")
	if "ideafix" in ffpe_filter:
		ffpe_filter_res_path = f"{source_set_name}/{dataset}/{ffpe_filter.lower().split("-")[0]}/{sample}/{sample}.{ffpe_filter.lower()}{res_suffix}"
	else:
		ffpe_filter_res_path = f"{source_set_name}/{dataset}/{ffpe_filter.lower()}/{sample}/{sample}.{ffpe_filter.lower()}{res_suffix}"


	ffpe_filter_res = pl.read_csv(ffpe_filter_res_path, separator="\t")
	ffpe_filter_res_micr_filtered = ffpe_filter_res.join(micr_artifacts, on =["chrom", "pos", "ref", "alt"], how="anti")

	output_dir = os.path.dirname(ffpe_filter_res_path).replace(source_set_name, new_set_name) # f"{vcf_filter_type.replace("filtered", "filtered_micr")}/{dataset}/{ffpe_filter.lower().split('-')[0]}/{sample}"
	os.makedirs(output_dir, exist_ok=True)
	output_path = f"{output_dir}/{sample}.{ffpe_filter.lower()}{res_suffix}"

	ffpe_filter_res_micr_filtered.write_csv(output_path, separator="\t")
	print(f"\t\tMICR filtered {ffpe_filter.upper()} results written to: {output_path}")

	variants_before = ffpe_filter_res.height
	artifact_variants = ffpe_filter_res.height - ffpe_filter_res_micr_filtered.height
	variants_after = ffpe_filter_res_micr_filtered.height
	print(f"\t\t{artifact_variants} MICR artifacts removed from {variants_before} variants ({(artifact_variants) / variants_before * 100:.2f}%).")
	print(f"\t\t{variants_after} variants remain after MICR filtering.")

	summary = pl.DataFrame({
		"source_sample_set": source_set_name,
		"sample": sample,
		"ffpe_filter" : f"{ffpe_filter.lower()}",
		"total_variants": variants_before,
		"micr_artifacts_removed": artifact_variants,
		"variants_after_micr_filtering": variants_after,
		"percent_removed": artifact_variants / variants_before * 100,
	})

	return summary



def filter_sample_set(
		msec_filter_col: str,
		source_set_name:str,
		new_set_name: str,
		dataset: str,
) -> None:
	
	print(f"Filtering {source_set_name}...")

	## Filter MICR artifacts from each dataset
	microsec_paths = glob.glob(f"{source_set_name}/{dataset}/microsec/*/*.microsec.tsv")
	# print(len(microsec_paths))
		
	micr_filtering_summaries = []

	for path in tqdm(microsec_paths):
		sample = os.path.basename(path).replace(".microsec.tsv", "")
		# print(f"\n{i+1}. Processing... Dataset: {source_set_name}, Sample: {sample}")
		
		microsec = pl.read_csv(path, separator="\t", infer_schema_length=1000).rename({"Chr": "chrom"}).rename(lambda x: x.lower())
		microsec_artifacts = microsec.filter(pl.col(msec_filter_col) == "Artifact suspicious")
		
		mobsnvf_summary = filter_micr_artifacts("mobsnvf", dataset, new_set_name, source_set_name, sample, microsec_artifacts)
		micr_filtering_summaries.append(mobsnvf_summary)
		
		sobdetector_summary = filter_micr_artifacts("sobdetector", dataset, new_set_name, source_set_name, sample, microsec_artifacts)
		micr_filtering_summaries.append(sobdetector_summary)

		vafsnvf_summary = filter_micr_artifacts("vafsnvf", dataset, new_set_name, source_set_name, sample, microsec_artifacts)
		micr_filtering_summaries.append(vafsnvf_summary)

		ideafix_xgb_summary = filter_micr_artifacts("ideafix-xgboost", dataset, new_set_name, source_set_name, sample, microsec_artifacts, res_suffix=".tsv")
		micr_filtering_summaries.append(ideafix_xgb_summary)

		gatk_obmm_summary = filter_micr_artifacts("gatk-obmm", dataset, new_set_name, source_set_name, sample, microsec_artifacts, res_suffix=".tsv")
		micr_filtering_summaries.append(gatk_obmm_summary)

		ffpolish_summary = filter_micr_artifacts("ffpolish", dataset, new_set_name, source_set_name, sample, microsec_artifacts, res_suffix=".tsv")
		micr_filtering_summaries.append(ffpolish_summary)
		
	micr_filtering_summaries_df = pl.concat(micr_filtering_summaries, how="diagonal_relaxed")
	micr_filtering_summaries_df.write_csv(f"{new_set_name}/{dataset}/micr_filtering_summary.tsv", separator="\t")
	

# %% [markdown]
# #### DP>=10

## Remove artifacts identified by MicroSEC filter 1,2,3,4
filter_sample_set(
	msec_filter_col = "msec_filter_1234", #"msec_filter_all"
    source_set_name = "mutect2-matched-normal_pass-orientation-dp-filtered",
	new_set_name = "mutect2-matched-normal_pass-orientation-dp-filtered_micr1234",
    dataset="FFX"
)

## Remove artifacts identified by all 8 filter MicroSEC
filter_sample_set(
	msec_filter_col = "msec_filter_all",
	source_set_name = "mutect2-matched-normal_pass-orientation-dp-filtered",
	new_set_name = "mutect2-matched-normal_pass-orientation-dp-filtered_micr",
	dataset="FFX"
)



# ## Remove artifacts identified by MicroSEC filter 1,2,3,4
# filter_sample_set(
# 	msec_filter_col = "msec_filter_1234", #"msec_filter_all"
#     source_set_name = "mutect2-matched-normal_pass-orientation-dp-filtered",
# 	new_set_name = "mutect2-matched-normal_pass-orientation-dp-filtered_micr1234",
#     dataset="FFG"
# )

# ## Remove artifacts identified by all 8 filter MicroSEC
# filter_sample_set(
# 	msec_filter_col = "msec_filter_all",
# 	source_set_name = "mutect2-matched-normal_pass-orientation-dp-filtered",
# 	new_set_name = "mutect2-matched-normal_pass-orientation-dp-filtered_micr",
# 	dataset="FFG"
# )


# %% [markdown]
# #### DP>=20

## Remove artifacts identified by MicroSEC filter 1,2,3,4
filter_sample_set(
	msec_filter_col = "msec_filter_1234", #"msec_filter_all"
    source_set_name = "mutect2-matched-normal_pass-orientation-dp20-filtered",
	new_set_name = "mutect2-matched-normal_pass-orientation-dp20-filtered_micr1234",
    dataset="FFX"
)

## Remove artifacts identified by all 8 filter MicroSEC
filter_sample_set(
	msec_filter_col = "msec_filter_all",
	source_set_name = "mutect2-matched-normal_pass-orientation-dp20-filtered",
	new_set_name = "mutect2-matched-normal_pass-orientation-dp20-filtered_micr",
	dataset="FFX"
)


# ## Remove artifacts identified by MicroSEC filter 1,2,3,4
# filter_sample_set(
# 	msec_filter_col = "msec_filter_1234", #"msec_filter_all"
#     source_set_name = "mutect2-matched-normal_pass-orientation-dp20-filtered",
# 	new_set_name = "mutect2-matched-normal_pass-orientation-dp20-filtered_micr1234",
#     dataset="FFG"
# )

# ## Remove artifacts identified by all 8 filter MicroSEC
# filter_sample_set(
# 	msec_filter_col = "msec_filter_all",
# 	source_set_name = "mutect2-matched-normal_pass-orientation-dp20-filtered",
# 	new_set_name = "mutect2-matched-normal_pass-orientation-dp20-filtered_micr",
# 	dataset="FFG"
# )



