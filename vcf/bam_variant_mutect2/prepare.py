#!/usr/bin/env python

# Prepare json input files for WDL mutect2 workflow 

import os
import json
import pandas as pd
import numpy as np
import glob
# ---


# absolute path is required in the json input files for WDL
ref_root = '../../data/reference_genome/GRCh38'
ref_fname = 'GRCh38.d1.vd1'
vcf_root = "../../data/gatk-best-practices/somatic-hg38"
bam_root = "../../data/bam"
bundle_root = "../../data/gatk-test-data/mutect2"
intervals_root = "../../data/ref"
gatk_path = "/home/moyukh/miniconda3/envs/ena-ffpe/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"

out_dir = 'inputs'

if not os.path.exists(out_dir):
	os.makedirs(out_dir)

wgs_normal_bam = os.path.abspath("../../data/bam/WGS/WGS_IL_N_1.bwa.dedup.bam")
wes_normal_bam = os.path.abspath("../../data/bam/WES/WES_IL_N_2.bwa.dedup.bam")
amplicon_normal_bam = os.path.abspath("../../data/bam/AmpliSeq_bams/AmpliSeq.bwa.HCC1395BL_1.bam")
wes_target_regions = os.path.abspath("../../data/reference_genome/Exome_Target_bed/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed")


ref_fpath = os.path.abspath(os.path.join(ref_root, ref_fname))
vcf_path = os.path.abspath(vcf_root)
bundle_path = os.path.abspath(bundle_root)
intervals_path = os.path.abspath(intervals_root)

# annot_path = '../../annot/sample_annotation.tsv'
# annot = pd.read_csv(annot_path, sep='\t')

bam_paths = glob.glob(f"{bam_root}/*/*.bam")

# base wdl input
base = {
	'bam_variant_mutect2.run_funcotator': False,
	'bam_variant_mutect2.ref_fasta': ref_fpath + '.fa',
	'bam_variant_mutect2.ref_fai': ref_fpath + '.fa.fai',
	'bam_variant_mutect2.ref_dict': ref_fpath + '.dict',
	'bam_variant_mutect2.pon': os.path.join(vcf_path, '1000g_pon.hg38.vcf.gz'),
	'bam_variant_mutect2.pon_idx': os.path.join(vcf_path, '1000g_pon.hg38.vcf.gz.tbi'),
	'bam_variant_mutect2.gnomad': os.path.join(vcf_path, 'af-only-gnomad.hg38.vcf.gz'),
	'bam_variant_mutect2.gnomad_idx': os.path.join(vcf_path, 'af-only-gnomad.hg38.vcf.gz.tbi'),
	'bam_variant_mutect2.variants_for_contamination':  os.path.join(vcf_path, 'small_exac_common_3.hg38.vcf.gz'),
	'bam_variant_mutect2.variants_for_contamination_idx': os.path.join(vcf_path, 'small_exac_common_3.hg38.vcf.gz.tbi'),
	'bam_variant_mutect2.realignment_index_bundle': os.path.join(bundle_path, 'Homo_sapiens_assembly38.index_bundle'),
	'bam_variant_mutect2.scatter_count': 4,
	'bam_variant_mutect2.gatk_docker': 'broadinstitute/gatk:4.6.2.0',
	'bam_variant_mutect2.gatk_override': gatk_path,
	'bam_variant_mutect2.bam_mutect2.mem': 4,
	'bam_variant_mutect2.run_orientation_bias_mixture_model_filter': True,
}


for path in bam_paths:
	
	study_alias = path.split("/")[-2]
	sample_name = os.path.basename(path).replace(".bam", "").replace(".dedup", "")
 
	if (("_N_" in sample_name) | ("HCC1395BL" in sample_name)):
		print(f"Skipping {sample_name} since it is normal")
		continue
	
	bam = path
	bai = f"{path}.bai"
	
	if not os.path.exists(bam):
		raise FileNotFoundError(f"{bam} does not exist")
	if not os.path.exists(bai):
		bai = f"{".".join(path.split(".")[0:-1])}.bai"
		if not os.path.exists(bai):
			raise FileNotFoundError(f"{bai} does not exist")

	# write wdl input json file for each sample
	out = base.copy()
	if study_alias in ["WES", "FFX"]:
		out['bam_variant_mutect2.intervals'] = os.path.join(intervals_path, 'wgs_calling_regions.hg38.interval_list')
		out['bam_variant_mutect2.normal_bam'] = wes_normal_bam
		out['bam_variant_mutect2.normal_bai'] = wes_normal_bam.replace(".bam", ".bai")
	elif study_alias in ["WGS", "FFG", "LBP"]:
		out['bam_variant_mutect2.intervals'] = os.path.join(intervals_path, 'wgs_calling_regions.hg38.interval_list')
		out['bam_variant_mutect2.normal_bam'] = wgs_normal_bam
		out['bam_variant_mutect2.normal_bai'] = wgs_normal_bam.replace(".bam", ".bai")
	elif study_alias in ["AmpliSeq_bams"]:
		out['bam_variant_mutect2.intervals'] = os.path.join(intervals_path, 'wgs_calling_regions.hg38.interval_list')
		out['bam_variant_mutect2.m2_extra_args'] = '--disable-read-filter NotDuplicateReadFilter --downsampling-stride 50 --linked-de-bruijn-graph --max-reads-per-alignment-start 0' #  --max-reads-per-alignment-start 500 --dont-use-soft-clipped-bases --annotations-to-exclude StrandBiasBySample --annotations-to-exclude ReadPosRankSumTest
		out['bam_variant_mutect2.normal_bam'] = amplicon_normal_bam
		out['bam_variant_mutect2.normal_bai'] = f"{amplicon_normal_bam}.bai"

	out['bam_variant_mutect2.tumor_bam'] = os.path.abspath(bam)
	out['bam_variant_mutect2.tumor_bai'] = os.path.abspath(bai)
	out_path = os.path.join(out_dir, f"{sample_name.replace(".", "-")}.inputs")

	with open(out_path, 'w') as outf:
		outf.write(json.dumps(out, indent=True, sort_keys=True))

	print(f"Prepared inputs for: {path}")

print("Done.")




