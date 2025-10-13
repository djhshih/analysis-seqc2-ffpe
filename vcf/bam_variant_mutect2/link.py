# %%
#!/usr/bin/env python3
import glob
import os
# import shutil

# %%
variant_paths = glob.glob("cromwell-executions/bam_variant_mutect2/*/call-vcf_filter/execution/*.vcf*")
outdir_root = "../mutect2-matched-normal_filtermutectcalls_obmm_unfiltered"

# %%
for i, path in enumerate(variant_paths):
    basename = os.path.basename(path)
    sample_id = basename.removesuffix(".idx").removesuffix("-filtered.vcf").split(".")[0]
    study = sample_id.split("_")[0]
    
    out_filename = basename.replace("-filtered", "")
    
    outdir = f"{outdir_root}/{study}/{sample_id}"
    os.makedirs(outdir, exist_ok=True)
    out_path = f"{outdir}/{out_filename}"
    
    if os.path.exists(out_path):
        print(f"'{out_filename}' already exists at '{outdir}'\n\tSKIPPING...")
        continue
    
    os.link(path, out_path)
    print(f"Created hard link for: {out_filename}")



