#!/usr/bin/env Rscript
library(io)
library(precrec)
library(dplyr)
source("../common-ffpe-snvf/R/eval.R")

## Setup
repo_root <- ".."
eval_dir = file.path(repo_root, "eval-wes-wgs-gt-union")
ffpe_snvf_dir = file.path(repo_root, "ffpe-snvf")

#######################################################

evaluate_sample_set <- function(
	model_name,
	dataset,
	variant_set,
	result_set = NULL,
	snvf_res_ext = "tsv"
) {

	if(is.null(result_set)){
		result_set <- variant_set
	}

	message(sprintf("Evaluating %s on Dataset: %s, Variant Set: %s", model_name, dataset, variant_set))

	if(grepl("ideafix", model_name)){
		snvf_paths <- Sys.glob(file.path(ffpe_snvf_dir, dataset, variant_set, unlist(strsplit(model_name, "-"))[1], sprintf("*/*.%s.%s", model_name, snvf_res_ext)))
	} else {
		snvf_paths <- Sys.glob(file.path(ffpe_snvf_dir, dataset, variant_set, model_name, sprintf("*/*.%s.%s", model_name, snvf_res_ext)))
	}
	snvf_paths <- Sys.glob(file.path(ffpe_snvf_dir, dataset, variant_set, unlist(strsplit(model_name, "-"))[1], sprintf("*/*.%s.%s", model_name, snvf_res_ext)))
	message(sprintf("	Found %d %s resuts...", length(snvf_paths), model_name))

	for (i in seq_along(snvf_paths)){

		snvf_path <- snvf_paths[i]
		sample_name <- basename(dirname(snvf_path))
		outdir_root <- file.path(dataset, result_set)

		message(sprintf("	%d. Sample: %s", i, sample_name))
		message(sprintf("		SNVF path: %s", snvf_path))

		## Read prepared ground truth
		ground_truth_path <- file.path("..", "ground-truth", dataset, gsub("-micr1234", "", variant_set), sample_name, sprintf("%s.ground-truth.tsv", sample_name))
		message(sprintf("		Reading Ground Truth from: %s", ground_truth_path))
		ground_truth <- read.delim(ground_truth_path)
		ground_truth$truth <- as.logical(ground_truth$truth)

		## Read model's score for current sample and apply model specific processing.
		## Then annotate ground truth
		d <- preprocess_filter(read.delim(snvf_path), model_name)
		d <- merge(ground_truth, d, by=c("chrom", "pos", "ref", "alt"))

		## Ensure that the ground truth annotation are not all true or all false
		## These if so they cannot be evaluated
		if((nrow(d[d$truth, ]) == 0)){
			message(sprintf("	no true labels exist for %s", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}
		if((nrow(d[!d$truth, ]) == 0)){
			message(sprintf("	no false labels exist for %s", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}

		# Evaluate the filter's performance
		res <- evaluate_filter(d, model_name, sample_name)

		# Write results
		write_sample_eval(d, res, outdir_root, sample_name, model_name)
		
	}

	# Overall Evaluation
	## The scores annotated with ground truth is combined into a single dataframe
	message("	Performing overall Evaluation across all samples")

	eval_data_paths <- Sys.glob(file.path(eval_dir, dataset, result_set, sprintf("model-scores_truths/*/*%s-scores_truths.tsv", model_name)))

	all_score_truth <- bind_rows(
		lapply(eval_data_paths, function(path) {
			d <- read.delim(path)
			d$sample_name <- sample_name
			return(d)
		})
	)

	# Evaluate across all samples
	overall_res <- evaluate_filter(all_score_truth, model_name, "all-samples")
	write_overall_eval(all_score_truth, overall_res, outdir_root, "all-samples", model_name)
	message(cat("\tComplete.\n"))

}

########################################################

model_name <- "ideafix-xgboost"

######################### SEQC2 FFX ###############################

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFX",
	variant_set = "mutect2-tn_filtered_pass-orientation"
)

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFX",
	variant_set = "mutect2-tn_filtered_pass-orientation-dp20"
)

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFX",
	variant_set = "mutect2-tn_filtered_pass-orientation-dp20-blacklist"
)

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFX",
	variant_set = "mutect2-tn_filtered_pass-orientation-exome"
)

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFX",
	variant_set = "mutect2-tn_filtered_pass-orientation-exome-dp20"
)

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFX",
	variant_set = "mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist"
)

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFX",
	variant_set = "mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist-clonal"
)


######################### SEQC2 FFG ###############################

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFG",
	variant_set = "mutect2-tn_filtered_pass-orientation"
)

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFG",
	variant_set = "mutect2-tn_filtered_pass-orientation-dp20"
)

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFG",
	variant_set = "mutect2-tn_filtered_pass-orientation-dp20-blacklist"
)

evaluate_sample_set(
	model_name = model_name,
	dataset = "FFG",
	variant_set = "mutect2-tn_filtered_pass-orientation-dp20-blacklist-clonal"
)

