#!/usr/bin/python3
import polars as pl
import pysam
import os
import glob


# %% [markdown]
# ### Function Definitions

## Create Mutation Info for MicroSEC

def import_formatted_vcf(vcf_path: str, ct_only: bool = False) -> pl.DataFrame:
	"""
	Reads a VCF file and returns a Polars DataFrame formatted for MicroSEC analysis.

	The function:
	- Reads a VCF file as a tab-delimited table, ignoring header lines starting with '##'.
	- Drops all contigs, only retaining Chr 1-22, X, and Y
	- Selects the relevant columns (#CHROM, POS, REF, ALT).
	- Adds placeholder columns required by MicroSEC: Sample, Mut_type, SimpleRepeat_TRF, and Neighborhood_sequence.
	- Renames columns to standardized names: Chr, Pos, Ref, Alt.
	- Splits multi-allelic ALT entries into separate rows, resulting in a biallelic representation for each variant.

	Parameters
	----------
	vcf_path : str
		Path to the VCF file to be imported.

	Returns
	-------
	pl.DataFrame
		A Polars DataFrame with columns: Sample, Mut_type, Chr, Pos, Ref, Alt, SimpleRepeat_TRF, Neighborhood_sequence.
		Each row represents a single biallelic variant.
	"""
	allowed_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
	vcf = (
		pl.read_csv(vcf_path, separator="\t", comment_prefix="##")
		.filter(pl.col("#CHROM").is_in(allowed_chroms))
		.select(["#CHROM", "POS", "REF", "ALT"])
		.with_columns([
			pl.lit("-").alias("Sample"),
			pl.lit("-").alias("Mut_type"),
			pl.lit("-").alias("SimpleRepeat_TRF"),
			pl.lit("-").alias("Neighborhood_sequence")
		])
		.rename({
			"#CHROM" : "Chr",
			"POS" : "Pos", 
			"REF": "Ref", 
			"ALT" : "Alt"
		})
		.select(["Sample", "Mut_type", "Chr", "Pos", "Ref", "Alt", "SimpleRepeat_TRF", "Neighborhood_sequence"])
  		.with_columns(
			pl.col("Alt").str.split(",")
		)
		.explode("Alt")
	)
 
	if(ct_only):
		vcf = vcf.filter(
			((pl.col("Ref") == "C") & (pl.col("Alt") == "T")) |
			((pl.col("Ref") == "G") & (pl.col("Alt") == "A"))
		)

	return vcf


def make_mut_info(vcf_path: str, sample_name: str, genome: pysam.FastaFile, ct_only: bool = False) -> pl.DataFrame:
    
	mut_info = import_formatted_vcf(vcf_path, ct_only).with_columns(pl.lit(sample_name).alias("Sample"))
	
	# Define a function to be applied to each row-like struct
	# Used to retrieve neighbourhood sequence via map_elements which is vectorized for speed
	def get_neighborhood_sequence(row_struct: dict) -> str:
		chrom = row_struct["Chr"]
		pos = row_struct["Pos"]
		alt = row_struct["Alt"]

		try:
			chr_length = genome.get_reference_length(chrom)

			# Calculate start and end for fetching flanks, ensuring they are within bounds
			left_start = max(0, pos - 21)
			left_end = max(0, pos - 1)

			right_start = pos
			right_end = min(chr_length, pos + 20)

			left_flank = genome.fetch(chrom, left_start, left_end)
			right_flank = genome.fetch(chrom, right_start, right_end)

			return (left_flank + alt + right_flank).upper()

		except ValueError:
			# Pysam can raise ValueError for contigs not in the FASTA file
			return None

	# Create a struct of the columns needed for the neighborhood sequence calculation
	# and then apply the function.
	mut_info = mut_info.with_columns(
		pl.struct(["Chr", "Pos", "Alt"])
		.map_elements(get_neighborhood_sequence, return_dtype=pl.Utf8)
		.alias("Neighborhood_sequence")
	)
	
	snv_mask = mut_info["Ref"].str.len_chars() == mut_info["Alt"].str.len_chars()
	deletion_mask = (mut_info["Ref"].str.len_chars() > 1) & (mut_info["Alt"].str.len_chars() == 1)
	insertion_mask = (mut_info["Ref"].str.len_chars() == 1) & (mut_info["Alt"].str.len_chars() > 1)
	
	mut_info = (
		mut_info.with_columns([
			pl.when(snv_mask)
			.then(mut_info["Ref"].str.len_chars().cast(pl.Utf8) + "-snv")
			.when(deletion_mask)
			.then((mut_info["Ref"].str.len_chars() - 1).cast(pl.Utf8) + "-del")
			.when(insertion_mask)
			.then((mut_info["Alt"].str.len_chars() - 1).cast(pl.Utf8) + "-ins")
			.otherwise(mut_info["Mut_type"])
			.alias("Mut_type")
		])
	)
	
	# MicroSEC requires if variant to be at least 200 nucleotides away from the terminal ends of the chromosome
	mut_info = mut_info.with_columns(
		pl.col("Chr").map_elements(lambda chr: genome.get_reference_length(chr), return_dtype=pl.Int64).alias("chr_length")
	).filter(
		(pl.col("Pos") > 200) & (pl.col("Pos") < (pl.col("chr_length") - 200))
	).drop("chr_length")
 
	return mut_info

# %% [markdown]
# ### FFX

## SEQC2 - FFX

print("Making inputs for SEQC2 - FFX")
ref_path = os.path.abspath("../../data/seqc2-reference-genome/GRCh38/GRCh38.d1.vd1.fa")
outdir_root = "../mutect2-matched-normal_pass-orientation-filtered/FFX/microsec/inputs"
genome = pysam.FastaFile(ref_path)


vcf_paths = sorted(glob.glob("../../vcf/mutect2-matched-normal_pass-orientation-filtered/FFX/*/*.vcf"))
bam_paths = sorted(glob.glob(f"../../data/bam/FFX/*.bam"))


# lookup_table_ffpe = lookup_table.filter(pl.col("preservation") == "FFPE")

vcf_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).split(".")[0] for path in vcf_paths],
	"vcf_path" : [os.path.abspath(path) for path in vcf_paths]
})

bam_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).split(".")[0] for path in bam_paths],
	"bam_path" : [os.path.abspath(path) for path in bam_paths]
})

vcf_bam_table = vcf_table.join(bam_table, on="sample_name")

for i in range(vcf_bam_table.shape[0]):

	sample_name = vcf_bam_table[i, "sample_name"]
	 
	outdir = f"{outdir_root}/mut_info"
	os.makedirs(outdir, exist_ok=True)

	mut_info = make_mut_info(vcf_bam_table[i, "vcf_path"], sample_name, genome, ct_only = False)

	mut_info.write_csv(f"{outdir}/{sample_name}.microsec.mut-info.tsv", separator="\t")
	print(f"Finished creating Mut-Info for - {i + 1}. {sample_name}")

### Create Sample Info for MicroSEC (polars)

mut_info_paths = sorted(glob.glob(f"{outdir_root}/mut_info/*.microsec.mut-info.tsv"))
mut_info_suffix = ".microsec.mut-info.tsv"

rows = []
for i, p in enumerate(mut_info_paths):
	sample_name = os.path.basename(p).split(mut_info_suffix)[0]

	bam_path = vcf_bam_table.filter(pl.col("sample_name") == sample_name).get_column("bam_path")[0]

	# These were the read lengths observed in the bam files
	if "_GZ_" in sample_name:
		read_length = 150
	elif "_IL_" in sample_name:
		read_length = 125

	rows.append({
		"sample_name": sample_name,
		"mutation information tsv file": os.path.abspath(p),
		"BAM file": bam_path,
		"read length": read_length,  # majority read length in bam files
		"adapter sequence read 1": None,
		"optional: adapter sequence read 2": None,
		"sample type: Human or Mouse": "hg38",
		"panel name": "TOP",
		"optional: reference genome fasta file": ref_path
	})

sample_info = pl.DataFrame(rows)

for i, sample_name in enumerate(sample_info.get_column("sample_name")):
    sample_info[i].write_csv(f"{outdir_root}/{sample_name}.microsec.sample_info.tsv", separator="\t", include_header=False)

sample_info_paths = glob.glob(f"{outdir_root}/*.microsec.sample_info.tsv")

script_dir = "script_microsec"
os.makedirs(script_dir, exist_ok=True)

for path in sample_info_paths:

	sample_name = os.path.basename(path).replace(".microsec.sample_info.tsv", "")
	script_path = f"{script_dir}/{sample_name}.microsec.sh"
	
	with open(script_path, "w") as file:
		file.writelines("#!/usr/bin/env bash\n")
		file.writelines(f"Rscript ../MicroSEC_serial.R ../{"/".join(path.split("/")[:4])} {"/".join(path.split("/")[4:])} N\n")

# %% [markdown]
# ### FFG

## SEQC2 - FFG

print("Making inputs for SEQC2 - FFG")
ref_path = os.path.abspath("../../data/seqc2-reference-genome/GRCh38/GRCh38.d1.vd1.fa")
outdir_root = "../mutect2-matched-normal_pass-orientation-filtered/FFG/microsec/inputs"
genome = pysam.FastaFile(ref_path)


vcf_paths = sorted(glob.glob("../../vcf/mutect2-matched-normal_pass-orientation-filtered/FFG/*/*.vcf"))
bam_paths = sorted(glob.glob(f"../../data/bam/FFG/*.bam"))


# lookup_table_ffpe = lookup_table.filter(pl.col("preservation") == "FFPE")

vcf_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).split(".")[0] for path in vcf_paths],
	"vcf_path" : [os.path.abspath(path) for path in vcf_paths]
})

bam_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).split(".")[0] for path in bam_paths],
	"bam_path" : [os.path.abspath(path) for path in bam_paths]
})

vcf_bam_table = vcf_table.join(bam_table, on="sample_name")

for i in range(vcf_bam_table.shape[0]):

	sample_name = vcf_bam_table[i, "sample_name"]
	 
	outdir = f"{outdir_root}/mut_info"
	os.makedirs(outdir, exist_ok=True)

	mut_info = make_mut_info(vcf_bam_table[i, "vcf_path"], sample_name, genome, ct_only = False)

	mut_info.write_csv(f"{outdir}/{sample_name}.microsec.mut-info.tsv", separator="\t")
	print(f"Finished creating Mut-Info for - {i + 1}. {sample_name}")

### Create Sample Info for MicroSEC

mut_info_paths = sorted(glob.glob(f"{outdir_root}/mut_info/*.microsec.mut-info.tsv"))
mut_info_suffix = ".microsec.mut-info.tsv"

rows = []
for i, p in enumerate(mut_info_paths):
	sample_name = os.path.basename(p).split(mut_info_suffix)[0]

	bam_path = vcf_bam_table.filter(pl.col("sample_name") == sample_name).get_column("bam_path")[0]

	# These were the read lengths observed in the bam files
	read_length = 151

	rows.append({
		"sample_name": sample_name,
		"mutation information tsv file": os.path.abspath(p),
		"BAM file": bam_path,
		"read length": read_length,  # majority read length in bam files
		"adapter sequence read 1": None,
		"optional: adapter sequence read 2": None,
		"sample type: Human or Mouse": "hg38",
		"panel name": "TOP",
		"optional: reference genome fasta file": ref_path
	})

sample_info = pl.DataFrame(rows)

for i, sample_name in enumerate(sample_info.get_column("sample_name")):
    sample_info[i].write_csv(f"{outdir_root}/{sample_name}.microsec.sample_info.tsv", separator="\t", include_header=False)



sample_info_paths = glob.glob(f"{outdir_root}/*.microsec.sample_info.tsv")

script_dir = "script_microsec"
os.makedirs(script_dir, exist_ok=True)

for path in sample_info_paths:

	sample_name = os.path.basename(path).replace(".microsec.sample_info.tsv", "")
	script_path = f"{script_dir}/{sample_name}.microsec.sh"

	with open(script_path, "w") as file:
		file.writelines("#!/usr/bin/env bash\n")
		file.writelines(f"Rscript ../MicroSEC_serial.R ../{"/".join(path.split("/")[:4])} {"/".join(path.split("/")[4:])} N\n")


