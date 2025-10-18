#!/usr/bin/env bash
library(io)
library(precrec)

source("../../common-ffpe-snvf/R/eval.R")

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
	search_path <- file.path(ffpe_snvf_dir, model_name, "*", sprintf("*.%s.snv", model_name))
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
		message(sprintf("	%s", sample_name))

		# Read and preprocess the data, annotating with truth labels
		d <- read.delim(path)
		d <- preprocess_mobsnvf(d, ground_truth_variants)

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

#################################################################################

model_name <- "mobsnvf"
message(sprintf("Evaluating %s: ", model_name))

################################# SEQC2 FFX  ########################################


# Setup Directories and lookup table for FFX dataset

dataset_id <- "FFX"
ff_dataset_id <- "WES"
message(sprintf("Dataset: %s", dataset_id))

# Directory for the Somatic VCFs
vcf_dir <- "../vcf/mutect2-matched-normal_pass-orientation-filtered"
# Directory for the FFPE SNV filtering results
ffpe_snvf_dir <- sprintf("../ffpe-snvf/mutect2-matched-normal_pass-orientation-filtered/%s", dataset_id)
# Root output directory
outdir_root <- sprintf("mutect2-matched-normal_pass-orientation-filtered/%s", dataset_id)
# Specific output directory for combined scores and truths
score_truth_outdir <- sprintf("%s/model-scores_truths", outdir_root)
# Specific output directory for evaluation plots and metrics
eval_outdir <- sprintf("%s/roc-prc-auc/precrec", outdir_root)


# Create a unified ground truth by taking the union of all variants from Fresh Frozen samples
message(sprintf("Generating unified ground truth from all Fresh Frozen samples (%s dataset)...", ff_dataset_id))
ff_paths <- Sys.glob(file.path(vcf_dir, sprintf("%s/*/*_T_*.vcf", ff_dataset_id)))
ff_variants <- snv_union(ff_paths)
message("Ground truth generated.")

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
    overall_res <- evaluate_filter(all_score_truth, model_name, "all-ffpe-wes-samples")

    ## Write overall result to disk
    write_overall_eval(all_score_truth, overall_res, score_truth_outdir, eval_outdir, "all-ffpe-wes-samples", model_name)
    message("Overall evaluation complete.")
} else {
    message("Skipping overall evaluation as no combined score/truth data was generated.")
}

################################# SEQC2 FFG  ########################################


# Setup Directories and lookup table for FFG dataset

dataset_id <- "FFG"
ff_dataset_id <- "WGS"
message(sprintf("Dataset: %s", dataset_id))

# Directory for the Somatic VCFs
vcf_dir <- "../vcf/mutect2-matched-normal_pass-orientation-filtered"
# Directory for the FFPE SNV filtering results
ffpe_snvf_dir <- sprintf("../ffpe-snvf/mutect2-matched-normal_pass-orientation-filtered/%s", dataset_id)
# Root output directory
outdir_root <- sprintf("mutect2-matched-normal_pass-orientation-filtered/%s", dataset_id)
# Specific output directory for combined scores and truths
score_truth_outdir <- sprintf("%s/model-scores_truths", outdir_root)
# Specific output directory for evaluation plots and metrics
eval_outdir <- sprintf("%s/roc-prc-auc/precrec", outdir_root)


# Create a unified ground truth by taking the union of all variants from Fresh Frozen samples
message(sprintf("Generating unified ground truth from all Fresh Frozen samples (%s dataset)...", ff_dataset_id))
ff_paths <- Sys.glob(file.path(vcf_dir, sprintf("%s/*/*_T_*.vcf", ff_dataset_id)))
ff_variants <- snv_union(ff_paths)
message("Ground truth generated.")

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
    overall_res <- evaluate_filter(all_score_truth, model_name, "all-ffpe-wgs-samples")

    ## Write overall result to disk
    write_overall_eval(all_score_truth, overall_res, score_truth_outdir, eval_outdir, "all-ffpe-wgs-samples", model_name)
    message("Overall evaluation complete.")
} else {
    message("Skipping overall evaluation as no combined score/truth data was generated.")
}
