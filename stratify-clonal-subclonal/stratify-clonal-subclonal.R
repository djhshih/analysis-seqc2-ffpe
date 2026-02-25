#!usr/bin/env Rscript
library(io)
source("../common-ffpe-snvf/R/macni_somatic.R")

repo_root <- ".."

#' Filter VCFs into clonal and subclonal mutation tables using MACNI somatic model.
#'
#' For each sample VCF in repo_root/vcf/<variant_set>/<dataset>, runs run_macni(),
#' classifies mutations as clonal or subclonal, and writes per-sample
#' <sample>.clonal.tsv and <sample>.subclonal.tsv files (chrom,pos,ref,alt).
#'
#' @param dataset Character. Dataset name.
#' @param variant_set Character. Variant set name.
#' @param outdir_root Character or NULL. Output root directory; if NULL a default
#'        path under repo_root/stratify-clonal-subclonal is used.
#' @param vcf_ext Character. VCF file extension (default "vcf").
#' @return Invisibly NULL. Side effect: creates output files.
filter_dataset <- function(dataset, variant_set, outdir_root=NULL, vcf_ext="vcf.gz") {

	message(sprintf("Stratifying dataset: %s, variant set: %s", dataset, variant_set))

	if (is.null(outdir_root)){
		outdir_root <- file.path(repo_root, "stratify-clonal-subclonal", dataset, variant_set)
	}
	message(sprintf("	Output directory: %s", outdir_root))

	vcf_paths <- Sys.glob(file.path(repo_root, "vcf", dataset, variant_set, sprintf("*/*.%s", vcf_ext)))
	message(sprintf("	Found %d VCF files to process.", length(vcf_paths)))

	for (i in seq_along(vcf_paths)){
		path <- vcf_paths[i]
		sample_name <- basename(dirname(path))
		message(sprintf("		%s. Processing sample: %s", i, sample_name))

		res <- run_macni(path, sample_name)
		res$is_clonal <- ifelse(with(res, (macni_pp > 0.5 & !(is.nan(macni_pp) & vaf == 1))), FALSE, TRUE)
		
		clonal_muts <- res[res$is_clonal, c("chrom", "pos", "ref", "alt")]
		subclonal_muts <- res[!res$is_clonal, c("chrom", "pos", "ref", "alt")]

		outdir <- file.path(outdir_root, sample_name)
		dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

		out_stem <- file.path(outdir, sample_name)
		write.table(clonal_muts, paste0(out_stem, ".clonal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
		write.table(subclonal_muts, paste0(out_stem, ".subclonal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
		
	}
	
	message(cat("\tComplete.\n\n"))
}

## WES FFPE Dataset
filter_dataset(
	dataset="FFX", 
	variant_set="mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist", 
)

## WES Frozen Dataset
filter_dataset(
	dataset="WES", 
	variant_set="mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist", 
)

## WGS FFPE Dataset
filter_dataset(
	dataset="FFG", 
	variant_set="mutect2-tn_filtered_pass-orientation-dp20-blacklist", 
)

## WGS Frozen Dataset
filter_dataset(
	dataset="WGS", 
	variant_set="mutect2-tn_filtered_pass-orientation-dp20-blacklist", 
)

