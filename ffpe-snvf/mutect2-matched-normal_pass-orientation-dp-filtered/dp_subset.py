import polars as pl
import glob
import os

paths = sorted(
	glob.glob("../mutect2-matched-normal_pass-orientation-filtered/*/*/*/*.snv") +
	glob.glob("../mutect2-matched-normal_pass-orientation-filtered/*/microsec/outputs/*/*.tsv") +
	glob.glob("../mutect2-matched-normal_pass-orientation-filtered/*/ideafix/*/*ideafix-*.tsv")
)

paths = list(filter(lambda path: "input.snv" not in path, paths))


filter_type = os.path.basename(os.getcwd())

for i, path in enumerate(paths):
	
	print(f"{i+1}. Processing file: {path}")
	
	if "microsec" in path:
		model = os.path.basename(path).split(".")[-2]
		dataset = path.split("/")[-5]
	else:
		dataset = path.split("/")[-4]
		model = os.path.basename(path).split(".")[-2]
  
	sample_name = os.path.basename(path).split(".")[0]

	vcf_path = glob.glob(f"../../vcf/{filter_type}/{dataset}/{sample_name}/{sample_name}*.vcf")[0]

	try:
		vcf_df = pl.read_csv(vcf_path, separator="\t", comment_prefix="##", null_values=".", columns=["#CHROM", "POS", "REF", "ALT"], infer_schema_length=1000).rename({"#CHROM": "CHROM"}).rename(lambda x: x.lower())
		vcf_df = vcf_df.with_columns(pl.col("alt").str.split(",")).explode("alt")
	except Exception as e:
		raise RuntimeError(f"Failed to read VCF file at {vcf_path}: {e}")

	snvf = pl.read_csv(path, separator="\t", infer_schema_length=1000)

	if model == "microsec":
		subset = snvf.join(vcf_df, left_on=["Chr", "Pos", "Ref", "Alt"], right_on=["chrom", "pos", "ref", "alt"], how="semi")
	else:
		subset = snvf.join(vcf_df, on=["chrom", "pos", "ref", "alt"], how="semi")
  
	if model == "microsec":
		output_dir = f"{dataset}/{model}/outputs/{sample_name}"
		os.makedirs(output_dir, exist_ok=True)
		output_path = f"{output_dir}/{sample_name}.{model}.tsv"
		subset.write_csv(output_path, separator="\t")

	else:
		output_dir = f"{dataset}/{model.split('-')[0]}/{sample_name}"
		os.makedirs(output_dir, exist_ok=True)
		if "ideafix" in model:
			output_path = f"{output_dir}/{sample_name}.{model}.tsv"
		else:
			output_path = f"{output_dir}/{sample_name}.{model}.snv"
		subset.write_csv(output_path, separator="\t")
  
	print(f"\tWritten subsetted data to {output_path}\n")

print("All files processed.")


