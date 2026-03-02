#!/usr/bin/env python3

# Prepare json input files for WDL mutect2 workflow 

import os
import json
import pandas as pd
import glob
import sys

repo_root = "../.."
## Local Dependency
sys.path.append(f"{repo_root}/common-ffpe-snvf/python")

from common import return_path_if_exists

# --- Setup ---
outdir = 'inputs'
os.makedirs(outdir, exist_ok=True)

exec_outdir = "batch_execution"
os.makedirs(exec_outdir)

gatk_path = "/home/moyukh/miniconda3/envs/ffpe-bench/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"
wdl_path = return_path_if_exists(f"{repo_root}/common-ffpe-snvf/wdl/bam_variant_mutect2.wdl", abs=True)

# absolute path is required in the json input files for WDL
## Reference Genome
ref_root = return_path_if_exists(f'{repo_root}/data/seqc2-reference-genome/GRCh38')
ref_fname = 'GRCh38.d1.vd1'
ref_fpath = os.path.join(ref_root, ref_fname)

## BAMs
bam_root = f"{repo_root}/data/bam"
wgs_normal_bam = return_path_if_exists(f"{bam_root}/WGS/WGS_IL_N_1.bwa.dedup.bam", abs=True)
wes_normal_bam = return_path_if_exists(f"{bam_root}/WES/WES_IL_N_2.bwa.dedup.bam", abs=True)
amplicon_normal_bam = return_path_if_exists(f"{bam_root}/AMS/AMS_AB_N_1.bwa.bam", abs=True)

## Target Regions
intervals_path = return_path_if_exists(f"{repo_root}/data/gatk-reference-genome", abs=True)
std_chr_list = return_path_if_exists(f"{repo_root}/data/misc/standard_chromosomes.list", abs=True)
# wes_target_regions = return_path_if_exists(f"{repo_root}/data/seqc2-reference-genome/Exome_Target_bed/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed", abs=True)

## Resources
resource_path = return_path_if_exists(f"{repo_root}/data/gatk-best-practices/somatic-hg38", abs=True)
# bundle_path = return_path_if_exists(f"{repo_root}/data/gatk-test-data/mutect2", abs=True)


## ----------------------------------------

bam_paths = glob.glob(f"{bam_root}/*/*_T_*.bam")

# base wdl input
base = {
	'bam_variant_mutect2.run_funcotator': False,
	'bam_variant_mutect2.ref_fasta': return_path_if_exists(f"{ref_fpath}.fa"),
	'bam_variant_mutect2.ref_fai':  return_path_if_exists(f"{ref_fpath}.fa.fai"),
	'bam_variant_mutect2.ref_dict':  return_path_if_exists(f"{ref_fpath}.dict"),
	'bam_variant_mutect2.pon': return_path_if_exists(os.path.join(resource_path, '1000g_pon.hg38.vcf.gz')),
	'bam_variant_mutect2.pon_idx': return_path_if_exists(os.path.join(resource_path, '1000g_pon.hg38.vcf.gz.tbi')),
	'bam_variant_mutect2.gnomad': return_path_if_exists(os.path.join(resource_path, 'af-only-gnomad.hg38.vcf.gz')),
	'bam_variant_mutect2.gnomad_idx': return_path_if_exists(os.path.join(resource_path, 'af-only-gnomad.hg38.vcf.gz.tbi')),
	'bam_variant_mutect2.variants_for_contamination':  return_path_if_exists(os.path.join(resource_path, 'small_exac_common_3.hg38.vcf.gz')),
	'bam_variant_mutect2.variants_for_contamination_idx': return_path_if_exists(os.path.join(resource_path, 'small_exac_common_3.hg38.vcf.gz.tbi')),
	# 'bam_variant_mutect2.realignment_index_bundle': return_path_if_exists(os.path.join(bundle_path, 'Homo_sapiens_assembly38.index_bundle')),
	'bam_variant_mutect2.scatter_count': 32,
	'bam_variant_mutect2.gatk_docker': 'broadinstitute/gatk:4.6.2.0',
	'bam_variant_mutect2.gatk_override': gatk_path,
	'bam_variant_mutect2.mutect2_scatter_mem_gb': 8,
	'bam_variant_mutect2.run_orientation_bias_mixture_model_filter': True,
}


for i, path in enumerate(bam_paths, start=1):

	print(f"{i}. Preparing inputs for: {path}")
	
	dataset = path.split("/")[-2]
	sample_name = os.path.basename(path).split(".")[0]
	print(f"\t Dataset: {dataset} | Sample Name: {sample_name}")
	
	bam = path
	bai = f"{path}.bai"
	
	if not os.path.exists(bai):
		bai = f"{".".join(path.split(".")[0:-1])}.bai"
		if not os.path.exists(bai):
			raise FileNotFoundError(f"{bai} does not exist")

	# write wdl input json file for each sample
	out = base.copy()
	if dataset in ["WES", "FFX"]:
		out['bam_variant_mutect2.intervals'] = os.path.join(intervals_path, 'wgs_calling_regions.hg38.interval_list')
		out['bam_variant_mutect2.normal_bam'] = wes_normal_bam
		out['bam_variant_mutect2.normal_bai'] = wes_normal_bam.replace(".bam", ".bai")
	elif dataset in ["WGS", "FFG", "LBP"]:
		out['bam_variant_mutect2.intervals'] = os.path.join(intervals_path, 'wgs_calling_regions.hg38.interval_list')
		out['bam_variant_mutect2.normal_bam'] = wgs_normal_bam
		out['bam_variant_mutect2.normal_bai'] = wgs_normal_bam.replace(".bam", ".bai")
	elif dataset in ["AMS"]:
		out['bam_variant_mutect2.intervals'] = std_chr_list
		out['bam_variant_mutect2.m2_extra_args'] = '--disable-read-filter NotDuplicateReadFilter --downsampling-stride 50 --linked-de-bruijn-graph --max-reads-per-alignment-start 0' #  --max-reads-per-alignment-start 500 --dont-use-soft-clipped-bases --annotations-to-exclude StrandBiasBySample --annotations-to-exclude ReadPosRankSumTest
		out['bam_variant_mutect2.normal_bam'] = amplicon_normal_bam
		out['bam_variant_mutect2.normal_bai'] = f"{amplicon_normal_bam}.bai"

	out['bam_variant_mutect2.tumor_bam'] = return_path_if_exists(bam, abs=True)
	out['bam_variant_mutect2.tumor_bai'] = return_path_if_exists(bai, abs=True)

	outpath = os.path.join(outdir, f"{sample_name.replace(".", "-")}.inputs")
	with open(outpath, 'w') as outf:
		outf.write(json.dumps(out, indent=True, sort_keys=True))

	exec_file=f"{exec_outdir}/{sample_name}_mutect2.sh"
	with open(exec_file, "w") as exec_file:
		exec_file.write(f"cromwell run {wdl_path} -i {return_path_if_exists(outpath, abs=True)}\n")

print("Done.")


