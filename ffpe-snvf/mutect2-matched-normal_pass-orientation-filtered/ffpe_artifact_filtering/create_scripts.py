import os
import glob
import polars as pl


ffpe_datasets = ["FFG", "FFX"]
models = ["mobsnvf", "vafsnvf", "sobdetector"]


## VCFs called using mutect2 with 1 matched normal
filtered_outdir_root = os.path.abspath("../../ffpe-snvf/mutect2-matched-normal_pass-orientation-filtered")

outdir = "script"
os.makedirs(outdir, exist_ok=True)

vcf_paths = glob.glob("../../vcf/mutect2-matched-normal_pass-orientation-filtered/*/*/*.vcf")

vcf_table = pl.DataFrame(
	{
		"sample_name" : [os.path.basename(path).split(".")[0] for path in vcf_paths],
		"vcf_path" : [os.path.abspath(path) for path in vcf_paths]
	}
)

bam_paths = glob.glob("../../data/bam/*/*.bam")

bam_table = pl.DataFrame(
	{
		"sample_name" : [os.path.basename(path).split(".")[0] for path in bam_paths],
		"bam_path" : [os.path.abspath(path) for path in bam_paths]
	}
)

bam_vcf_table = bam_table.join(vcf_table, on="sample_name", how="inner")

for i in range(bam_vcf_table.shape[0]):
	
	bam_path = bam_vcf_table[i, "bam_path"]
	vcf_path = bam_vcf_table[i, "vcf_path"]
	sample_name = os.path.basename(vcf_path).replace(".vcf", "").replace(".gz", "")
	dset = bam_path.split("/")[-2]
 
	if dset not in ffpe_datasets:
		print(f"{sample_name} is not FFPE")
		continue
	
	for model in models:
    
		filtered_outdir = f"{filtered_outdir_root}/{dset}"
		
		content = f"#!/bin/bash\nbash ../{model}.sh.template '{bam_path}' '{vcf_path}' '{filtered_outdir}'\n"
		filename = f"{outdir}/{model}_{sample_name}.sh"
		
		with open(filename, "w") as f:
			f.write(content)
			
		print(f"Created {model} filtering script for {vcf_path}")



# ## Datasets from the SEQCII NCBI repository

# outdir = "script_seqc2_ftp"
# os.makedirs(outdir, exist_ok=True)

# vcf_paths = glob.glob("../../vcf/seqc2_ftp/*/*.vcf.gz")
# filtered_outdir_root = os.path.abspath("../../ffpe-snvf/seqc2_ftp")
# bam_dir = "../../data/bam"
# bam_paths = glob.glob(f"{bam_dir}/*/*.bam")

# # Filter by FFPE Datasets
# vcf_paths = [path for path in vcf_paths if any(dset in path for dset in ffpe_datasets)]

# vcf_table = pl.DataFrame({
# 	"sample_name": [os.path.basename(path).split(".")[0] for path in vcf_paths],
# 	"vcf_path": [os.path.abspath(path) for path in vcf_paths],
# })

# vcf_table = vcf_table.filter(
# 	pl.col("sample_name").str.contains_any(ffpe_datasets)
# ).with_columns(
# 	pl.col("sample_name").str.split("_VS_").list.get(0).alias("sample_name")
# ).sort("sample_name")

# bam_table = pl.DataFrame({
# 	"sample_name": [os.path.basename(path).split(".")[0] for path in bam_paths],
# 	"bam_path": [os.path.abspath(path) for path in bam_paths],
# })

# bam_table = bam_table.filter(
# 	pl.col("sample_name").str.contains_any(ffpe_datasets)
# ).sort("sample_name")

# bam_vcf_table = vcf_table.join(bam_table, on="sample_name")
# bam_vcf_table


# for i in range(bam_vcf_table.shape[0]):
	
# 	bam_path = bam_vcf_table[i, "bam_path"]
# 	vcf_path = bam_vcf_table[i, "vcf_path"]
# 	sample_name = os.path.basename(vcf_path).replace(".vcf", "").replace(".gz", "")
# 	dset = bam_path.split("/")[-2]
	
# 	for model in models:
    
# 		filtered_outdir = f"{filtered_outdir_root}/{dset}"
		
# 		content = f"#!/bin/bash\nbash ../{model}.sh.template '{bam_path}' '{vcf_path}' '{filtered_outdir}'\n"
# 		filename = f"{outdir}/{model}_{sample_name}.sh"
		
# 		with open(filename, "w") as f:
# 			f.write(content)
			
# 		print(f"Created {model} filtering script for {vcf_path}")