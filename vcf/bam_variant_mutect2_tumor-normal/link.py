#!/usr/bin/env python
import glob
import os

## Path setup
repo_root = "../.."
vcf_root = f"{repo_root}/vcf"
var_call_dir_name = "bam_variant_mutect2_tumor-normal"
variant_paths = glob.glob(f"{vcf_root}/{var_call_dir_name}/cromwell-executions/bam_variant_mutect2/*/call-vcf_filter/execution/*.vcf*")
out_variant_set = "mutect2-tumor-normal_filtermutectcalls_obmm_unfiltered"

## Link VCF and index
for i, path in enumerate(variant_paths, start=1):
    sample_name = os.path.basename(path).split(".")[0]
    dataset = sample_name.split("_")[0]
    
    out_filename = os.path.basename(path).replace("-filtered", "")
    
    outdir = f"{vcf_root}/{dataset}/{out_variant_set}/{sample_name}"
    os.makedirs(outdir, exist_ok=True)
    out_path = f"{outdir}/{out_filename}"
    
    if os.path.exists(out_path):
        print(f"'{out_filename}' already exists at '{outdir}'\n\tSKIPPING...")
        continue
    
    os.link(path, out_path)
    print(f"{i}. Created hard link for: {out_filename}")
