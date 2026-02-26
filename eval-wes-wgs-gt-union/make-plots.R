## This script is used to make ROC and PRC plots from the evaluation results aggregated using combine_results.R
source("../common-ffpe-snvf/R/plot.R")

dset_dirs <- c(
	"FFX/mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist-clonal",
	"FFX/mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist",
	"FFX/mutect2-tn_filtered_pass-orientation-exome-dp20",
	"FFX/mutect2-tn_filtered_pass-orientation-exome",
	"FFX/mutect2-tn_filtered_pass-orientation-dp20-blacklist",
	"FFX/mutect2-tn_filtered_pass-orientation-dp20",
	"FFX/mutect2-tn_filtered_pass-orientation"
	"FFG/mutect2-tn_filtered_pass-orientation-dp20-blacklist-clonal",
	"FFG/mutect2-tn_filtered_pass-orientation-dp20-blacklist",
	"FFG/mutect2-tn_filtered_pass-orientation-dp20",
	"FFG/mutect2-tn_filtered_pass-orientation"
)

dset_author <- c(
	"FFX" = "SEQC-II WES FFPE",
	"FFG" =  "SEQC-II WGS FFPE"
)

models <- c(
	"mobsnvf" = "MOBSNVF",
	"vafsnvf" = "VAFSNVF",
	"gatk-obmm" = "GATK-OBMM",
	"sobdetector" = "SOBDetector",
	"microsec" = "MicroSEC",
	"ideafix-xgboost" = "Ideafix",
	"ffpolish" = "FFPolish"
)

for (dir in dset_dirs){

	message(sprintf("Processing: %s", dir))

	## Set output directory
	outdir_root <- file.path(dir, "plots")

	## Make a vector with paths to all the precrec eval objects
	roc_coord_paths <- sort(list.files(file.path(dir, "roc-prc-auc/precrec"), pattern = "all-models_roc_coordinates.tsv", recursive = TRUE, full.names = TRUE))
	prc_coord_paths <- sort(list.files(file.path(dir, "roc-prc-auc/precrec"), pattern = "all-models_prc_coordinates.tsv", recursive = TRUE, full.names = TRUE))

	## Create plots for each of the evaluated samples
	message("Creating ROC PRC plot for:")

	outdir <- glue("{outdir_root}/roc_prc_plots")
	dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

	for (i in seq_len(length(roc_coord_paths))) {

		sample_name <- match_return_sample_name(roc_coord_paths[i], prc_coord_paths[i])
		message(glue("\t {i}. {sample_name}"))

		roc_coord <- qread(roc_coord_paths[i])
		roc_coord$model <- ifelse(roc_coord$model %in% names(models), models[roc_coord$model], roc_coord$model)

		prc_coord <- qread(prc_coord_paths[i])
		prc_coord$model <- ifelse(prc_coord$model %in% names(models), models[prc_coord$model], prc_coord$model)

		plots <- make_roc_prc_plot(roc_coord, prc_coord, title = NULL, subtitle = NULL, caption = NULL, text_scale=2, line_width = 1.5, legend_rows = 3, individual_plots = TRUE, legend_scale = 0.8)

		qdraw(plots$roc_prc, glue("{outdir}/{sample_name}_roc_prc_plot.pdf"), width = 8, height = 6)
		dir.create(glue("{outdir}/roc"), recursive = TRUE, showWarnings = FALSE)
		qdraw(plots$roc, glue("{outdir}/roc/{sample_name}_roc_plot.pdf"), width = 5, height = 6)
		dir.create(glue("{outdir}/prc"), recursive = TRUE, showWarnings = FALSE)
		qdraw(plots$prc, glue("{outdir}/prc/{sample_name}_prc_plot.pdf"), width = 5, height = 6)

	}
}