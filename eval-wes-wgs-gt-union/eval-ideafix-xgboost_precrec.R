#!/usr/bin/env bash
library(io)
library(precrec)

source("../common-ffpe-snvf/R/eval.R")

#####################################################################################

# Function for performing per-sample evaluation using a unified ground truth derived
# from the union of all Fresh Frozen (FF) samples.
# @param ffpe_snvf_dir (string) | Directory containing the model's SNV filtering results for a specific dataset.
# @param model_name (string) | Name of the model being evaluated, used in file paths.
# @param ground_truth_variants (data.frame) | A data frame of all true variants, typically created by snv_union().
# @param outdir_root (string) | Root output directory for writing evaluation results.

process_samples <- function(ffpe_snvf_dir, model_name, ground_truth_variants, outdir_root) {
	message("Processing samples with a unified ground truth:")

	# Construct the search path to find all of the model's result files
	search_path <- file.path(ffpe_snvf_dir, unlist(strsplit(model_name, "-")), "*", sprintf("*.%s.tsv", model_name))
	model_result_paths <- Sys.glob(search_path)

	if (length(model_result_paths) == 0) {
		message(sprintf("Warning: No result files found for model '%s' at search path: %s", model_name, search_path))
		return()
	}

	message(sprintf("Found %d result files to process.", length(model_result_paths)))

	# Evaluate each sample against the unified ground truth
	for (path in model_result_paths) {
		# Extract sample name from the file path. Assumes format 'sample_name.model_name.snv'
		sample_name <- unlist(strsplit(basename(path), "\\."))[1]
		message(sprintf("	%s: %s", sample_name, path))

		# Read and preprocess the data, annotating with truth labels
		d <- read.delim(path)
		d <- preprocess_ideafix(d, ground_truth_variants)

		# Check if truth labels are not exclusively TRUE or FALSE.
		# Cases like these are skipped as evaluation is not supported by precrec.
		if (nrow(d[d$truth, ]) == 0 | nrow(d[!d$truth, ]) == 0) {
			message(sprintf("		truth labels consists of only one class for %s. Skipping evaluation for this sample.", sample_name))
			next
		}

		# Evaluate the filter's performance for the individual sample
		res <- evaluate_filter(d, model_name, sample_name)

		# Write per-sample results (scores/truths table and evaluation metrics)
		write_sample_eval(d, res, outdir_root, sample_name, model_name)
	}
}


# Function for combining the ground truth annotated SNV score tables from a directory structure.
# @param score_truth_outdir (string) | Directory where the per-sample ground truth annotated scores were saved.
# @param model_name (string) | Name of model being evaluated. Used to find the correct files.

combine_snv_score_truth <- function(score_truth_outdir, model_name) {
	message("Combining all per-sample ground truth SNV score tables into one")

	# Construct the search path for the score/truth files written by write_sample_eval
	search_pattern <- sprintf("*_%s-scores_truths.tsv", model_name)
	score_truth_paths <- Sys.glob(file.path(score_truth_outdir, "*", search_pattern))

	if (length(score_truth_paths) == 0) {
		message(sprintf("Warning: No score/truth files found in '%s' for model '%s'", score_truth_outdir, model_name))
		return(data.frame()) # Return empty data frame if no files are found
	}
    
	message(sprintf("Found %d score/truth files to combine.", length(score_truth_paths)))

	# Read each file, add the sample name, and combine into a single data frame
	all_score_truth <- do.call(
		rbind,
		lapply(score_truth_paths, function(path) {
			# Extract sample name from the file name
			file_basename <- basename(path)
			sample_name <- gsub(sprintf("_%s-scores_truths.tsv", model_name), "", file_basename)
			
			message(sprintf("	Reading data for sample %s", sample_name))
			
			df <- read.delim(path)
			df$sample_name <- sample_name
			df
		})
	)
	
	return(all_score_truth)
}

evaluate_dataset <- function(
	dataset_id,
	ff_dataset_id,
	ff_vcf_dir,
	ffpe_snvf_dir,
	outdir_root,
	model_name,
	agg_name
){

	message(sprintf("Dataset: %s", dataset_id))

	# Directory for the FFPE SNV filtering results
	ffpe_snvf_dir <- file.path(ffpe_snvf_dir, dataset_id)
	# Root output directory
	outdir_root <- file.path(outdir_root, dataset_id)
	# Specific output directory for combined scores and truths
	score_truth_outdir <- file.path(outdir_root, "model-scores_truths")
	# Specific output directory for evaluation plots and metrics
	eval_outdir <- file.path(outdir_root, "roc-prc-auc/precrec")

	##----------------

	# Create a unified ground truth by taking the union of all variants from Fresh Frozen samples
	message("Generating unified ground truth from all Fresh Frozen samples (WES + WGS dataset)...")
	ff_paths <- c(Sys.glob(file.path(ff_vcf_dir, "WGS/*/*_T_*.vcf")), Sys.glob(file.path(ff_vcf_dir, "WES/*/*_T_*.vcf")))

	dir.create(file.path("ground-truth", basename(ff_vcf_dir)), recursive = TRUE, showWarnings = FALSE)
	gt_path <- file.path("ground-truth", basename(ff_vcf_dir) ,"WES_WGS_ff-variants_set.tsv")
	if (file.exists(gt_path)){
		message("Using pregenerated ground truth")
		ff_variants <- qread(gt_path)
	} else {
		ff_variants <- snv_union(ff_paths)
		message("Ground truth generated.")
		qwrite(ff_variants, gt_path)
	}

	##------------------

	# Perform per-sample evaluation
	process_samples(
		ffpe_snvf_dir = ffpe_snvf_dir, 
		model_name = model_name, 
		ground_truth_variants = ff_variants, 
		outdir_root = outdir_root
	)

	# Combine the per-sample score/truth tables into a single data frame
	all_score_truth <- combine_snv_score_truth(
		score_truth_outdir = score_truth_outdir, 
		model_name = model_name
	)

	# Proceed with overall evaluation if combined data is available
	if (nrow(all_score_truth) > 0) {
		message("Performing overall evaluation...")
		## Evaluate overall result
		overall_res <- evaluate_filter(all_score_truth, model_name, agg_name)

		## Write overall result to disk
		write_overall_eval(all_score_truth, overall_res, outdir_root, agg_name, model_name)
		message("Overall evaluation complete.")
	} else {
		message("Skipping overall evaluation as no combined score/truth data was generated.")
	}
	
}

#################################################################################

model_name <- "ideafix-xgboost"
message(sprintf("Evaluating %s: ", model_name))

# ################################# SEQC2 FFX w/ Matched Normal ########################################

# # Filtered FFX dataset vs Unfiltered WES VCF 
# evaluate_dataset(
# 	dataset_id = "FFX",
# 	ff_dataset_id = "WES",
# 	ff_vcf_dir = "../vcf/mutect2-matched-normal_filtermutectcalls_obmm_unfiltered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-filtered",
# 	outdir_root = "mutect2-matched-normal_pass-orientation-filtered.vs.unfiltered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wes-samples"
# )


# # Filtered FFX dataset vs Filtered WES VCF 
# evaluate_dataset(
# 	dataset_id = "FFX",
# 	ff_dataset_id = "WES",
# 	ff_vcf_dir = "../vcf/mutect2-matched-normal_pass-orientation-filtered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-filtered",
# 	outdir_root = "mutect2-matched-normal_pass-orientation-filtered.vs.filtered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wes-samples"
# )

# # Filtered FFX dataset vs Filtered WES VCF, DP >= 10
# evaluate_dataset(
# 	dataset_id = "FFX",
# 	ff_dataset_id = "WES",
# 	ff_vcf_dir = "../vcf/mutect2-matched-normal_pass-orientation-dp-filtered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-dp-filtered",
# 	outdir_root = "mutect2-matched-normal_pass-orientation-dp-filtered.vs.filtered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wes-samples"
# )

# # Filtered FFX dataset vs Filtered WES VCF; DP>=20
# evaluate_dataset(
# 	dataset_id = "FFX",
# 	ff_dataset_id = "WES",
# 	ff_vcf_dir = "../vcf/mutect2-matched-normal_pass-orientation-dp20-filtered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-dp20-filtered",
# 	outdir_root = "mutect2-matched-normal_pass-orientation-dp20-filtered.vs.filtered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wes-samples"
# )

# # Filtered FFX dataset vs Filtered WES VCF; MICR removed
# evaluate_dataset(
# 	dataset_id = "FFX",
# 	ff_dataset_id = "WES",
# 	ff_vcf_dir = "../vcf/mutect2-matched-normal_pass-orientation-filtered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-filtered_micr",
# 	outdir_root = "mutect2-matched-normal_pass-orientation-filtered_micr.vs.filtered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wes-samples"
# )

# Filtered FFX dataset vs Filtered WES VCF, DP >= 10; MICR1234 removed
evaluate_dataset(
	dataset_id = "FFX",
	ff_dataset_id = "WES",
	ff_vcf_dir = "../vcf/mutect2-matched-normal_pass-orientation-dp-filtered",
	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-dp-filtered_micr1234",
	outdir_root = "mutect2-matched-normal_pass-orientation-dp-filtered_micr1234.vs.filtered-ff",
	model_name = model_name,
	agg_name = "all-ffpe-wes-samples"
)

# Filtered FFX dataset vs Filtered WES VCF; DP>=20; MICR1234 removed
evaluate_dataset(
	dataset_id = "FFX",
	ff_dataset_id = "WES",
	ff_vcf_dir = "../vcf/mutect2-matched-normal_pass-orientation-dp20-filtered",
	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-dp20-filtered_micr1234",
	outdir_root = "mutect2-matched-normal_pass-orientation-dp20-filtered_micr1234.vs.filtered-ff",
	model_name = model_name,
	agg_name = "all-ffpe-wes-samples"
)

# ################################# SEQC2 FFG w/ Matched Normal ########################################

# # FFG dataset vs Unfiltered WGS VCF 
# evaluate_dataset(
# 	dataset_id = "FFG",
# 	ff_dataset_id = "WGS",
# 	ff_vcf_dir = "../vcf/mutect2-matched-normal_filtermutectcalls_obmm_unfiltered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-filtered",
# 	outdir_root = "mutect2-matched-normal_pass-orientation-filtered.vs.unfiltered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wgs-samples"
# )

# # FFG dataset vs Filtered WGS VCF 
# evaluate_dataset(
# 	dataset_id = "FFG",
# 	ff_dataset_id = "WGS",
# 	ff_vcf_dir = "../vcf/mutect2-matched-normal_pass-orientation-filtered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-filtered",
# 	outdir_root = "mutect2-matched-normal_pass-orientation-filtered.vs.filtered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wgs-samples"
# )


# FFG dataset vs Filtered WGS VCF, DP>=10
evaluate_dataset(
	dataset_id = "FFG",
	ff_dataset_id = "WGS",
	ff_vcf_dir = "../vcf/mutect2-matched-normal_pass-orientation-dp-filtered",
	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-dp-filtered",
	outdir_root = "mutect2-matched-normal_pass-orientation-dp-filtered.vs.filtered-ff",
	model_name = model_name,
	agg_name = "all-ffpe-wgs-samples"
)

# FFG dataset vs Filtered WGS VCF, DP>=20
evaluate_dataset(
	dataset_id = "FFG",
	ff_dataset_id = "WGS",
	ff_vcf_dir = "../vcf/mutect2-matched-normal_pass-orientation-dp20-filtered",
	ffpe_snvf_dir = "../ffpe-snvf/mutect2-matched-normal_pass-orientation-dp20-filtered",
	outdir_root = "mutect2-matched-normal_pass-orientation-dp20-filtered.vs.filtered-ff",
	model_name = model_name,
	agg_name = "all-ffpe-wgs-samples"
)


# ################################# SEQC2 FFX Tumor Only ########################################
# message("Tumor-Only VCF evaluation")


# # Filtered FFX dataset vs Filtered WES VCF 
# evaluate_dataset(
# 	dataset_id = "FFX",
# 	ff_dataset_id = "WES",
# 	ff_vcf_dir = "../vcf/mutect2-tumor-only_pass-orientation-filtered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-tumor-only_pass-orientation-filtered",
# 	outdir_root = "mutect2-tumor-only_pass-orientation-filtered.vs.filtered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wes-samples"
# )


# # Filtered FFX dataset vs Unfiltered WES VCF 
# evaluate_dataset(
# 	dataset_id = "FFX",
# 	ff_dataset_id = "WES",
# 	ff_vcf_dir = "../vcf/mutect2-tumor-only_filtermutectcalls_obmm_unfiltered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-tumor-only_pass-orientation-filtered",
# 	outdir_root = "mutect2-tumor-only_pass-orientation-filtered.vs.unfiltered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wes-samples"
# )

################################# SEQC2 FFG Tumor Only ########################################

# # FFG dataset vs Unfiltered WGS VCF 
# evaluate_dataset(
# 	dataset_id = "FFG",
# 	ff_dataset_id = "WGS",
# 	ff_vcf_dir = "../vcf/mutect2-tumor-only_filtermutectcalls_obmm_unfiltered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-tumor-only_pass-orientation-filtered",
# 	outdir_root = "mutect2-tumor-only_pass-orientation-filtered.vs.unfiltered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wgs-samples"
# )

# # FFG dataset vs Filtered WGS VCF 
# evaluate_dataset(
# 	dataset_id = "FFG",
# 	ff_dataset_id = "WGS",
# 	ff_vcf_dir = "../vcf/mutect2-tumor-only_pass-orientation-filtered",
# 	ffpe_snvf_dir = "../ffpe-snvf/mutect2-tumor-only_pass-orientation-filtered",
# 	outdir_root = "mutect2-tumor-only_pass-orientation-filtered.vs.filtered-ff",
# 	model_name = model_name,
# 	agg_name = "all-ffpe-wgs-samples"
# )

