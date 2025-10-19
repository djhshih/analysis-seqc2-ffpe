## This script is used to make ROC and PRC plots from the evaluation results aggregated using combine_results.R
source("../common-ffpe-snvf/R/plot.R")

dset_dirs <- c(
	"mutect2-matched-normal_pass-orientation-filtered/FFX",
	"mutect2-matched-normal_pass-orientation-filtered/FFG"
)

dset_author <- c(
	"FFX" = "SEQC-II WES FFPE",
	"FFG" =  "SEQC-II WGS FFPE"
)

for (dir in dset_dirs){

	dset <- unlist(strsplit(dir, "/"))[2]
	message(sprintf("Processing dataset: %s", dset))

	## Set output directory
	outdir_root <- file.path(dir, "plots_publication")

	## Make a vector with paths to all the precrec eval objects
	roc_coord_paths <- sort(list.files(file.path(dir, "roc-prc-auc/precrec"), pattern = "all-models_roc_coordinates.tsv", recursive = TRUE, full.names = TRUE))
	prc_coord_paths <- sort(list.files(file.path(dir, "roc-prc-auc/precrec"), pattern = "all-models_prc_coordinates.tsv", recursive = TRUE, full.names = TRUE))

	## Plot display size parameters for debugging
	# options(repr.plot.width = 6, repr.plot.height = 5)

	## Create plots for each of the evaluated samples
	message("Creating ROC PRC plot for:")

	outdir <- glue("{outdir_root}/roc_prc_plots")
	dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

	for (i in seq_len(length(roc_coord_paths))) {

		sample_name <- match_return_sample_name(roc_coord_paths[i], prc_coord_paths[i])
		message(glue("\t {i}. {sample_name}"))

		# ## Read in the variant set to count number and annotate the number of SNVs in the plot
		# eval_snv_set <- qread(glue("PRJEB8754/vcf_pass-orient-pos-sb_ad_filtered/model-scores_truths/{sample_name}/{sample_name}_model-scores_truths.tsv"))
		# snv_count <- nrow(eval_snv_set)

		roc_coord <- qread(roc_coord_paths[i])
		prc_coord <- qread(prc_coord_paths[i])

		plots <- make_roc_prc_plot(roc_coord, prc_coord, title = dset_author[dset], subtitle = NULL, caption = NULL, text_scale=1.5, line_width = 1, individual_plots = TRUE)

		qdraw(plots$roc_prc, glue("{outdir}/{sample_name}_roc_prc_plot.pdf"), width = 7, height = 5)
		dir.create(glue("{outdir}/roc"))
		qdraw(plots$roc, glue("{outdir}/roc/{sample_name}_roc_plot.pdf"), width = 4, height = 5)
		dir.create(glue("{outdir}/prc"))
		qdraw(plots$prc, glue("{outdir}/prc/{sample_name}_prc_plot.pdf"), width = 4, height = 5)

	}
}