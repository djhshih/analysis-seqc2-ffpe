#!/usr/bin/env python
import os
import polars as pl
import sys

repo_root = "../.."
# Local dependencies
sys.path.append(f"{repo_root}/common-ffpe-snvf/python")
from common import read_variants


## Main
snv_gt_path = f"{repo_root}/data/release/latest/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
gt_outdir = f"{repo_root}/ground-truth/from-seqc"
os.makedirs(gt_outdir, exist_ok=True)

snv_gt = read_variants(snv_gt_path, columns = ["#CHROM", "POS", "REF", "ALT", "FILTER"])
snv_gt_pass_hc = snv_gt.filter(pl.col("filter") == "PASS;HighConf")

snv_gt.write_csv(f"{gt_outdir}/seqc2_ground-truth_snv_pass-high-med-conf.tsv", separator="\t")
snv_gt_pass_hc.write_csv(f"{gt_outdir}/seqc2_ground-truth_snv_pass-high-conf.tsv", separator="\t")
