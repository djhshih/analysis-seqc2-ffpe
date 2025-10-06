import polars as pl


meta = pl.read_csv("sample-info_stage0.tsv", separator="\t")

# %% [markdown]
# ## Information from Paper
# 
# We generated WGS and WES data using various NGS library preparation protocols, seven NGS platforms (NovaSeq, HiSeq, PacBio, 10X Genomics, Ion Torrent, Miseq, and Affymetrix CytoScan HD) at six centers including Illumina (IL), National Cancer Institute (NC), Novartis (NV), European Infrastructure for Translational Medicine (EA), Fudan University (FD), and Loma Linda University (LL) 
# 
# 

sequencing_center = {
    "EA": "EATRIS",
    "FD": "Fudan University",
    "GT": "Genentech",
    "GV": "Garvan",
    "HA": "HudsonAlpha",
    "IL": "Illumina",
    "LL": "Linda Loma University",
    "MC": "Mayo Clinic",
    "NC": "NCI",
    "NV": "Novartis",
    "WC": "Weil Cornell"
}

study_name = {
	"WGS": "WGS of Fresh DNA extracted from HCC1395BL and HCC1395 cell lines",
	"FFG": "WGS of FFPE with fixation time of 1, 2, 6, or 24 hours from HCC1395BL and HCC1395 cell lines",
	"WES": "WES of Fresh DNA extracted from HCC1395BL and HCC1395 cell lines",
	"FFX": "WES of FFPE with fixation time of 1, 2, 6, or 24 hours from HCC1395BL and HCC1395 cell lines",
	"LBP": "WGS Libraries were made from different library protocols such as TruSeq Nano, TruSeq PCR Free and Nextera Flex library protocol with different input amount and sequenced on Illumina HiSeq 3000/4000",
	"SPP": "WGS Libraries were made from pooling the HCC1395 and HCC1395BL cell lines with various ratios (3:1, 1:1, 1:4, 1:9 and 1:19) to create mixtures.",
	"CHR": "10X Genomics Chromium Genome Sequencing (10X WGS) data sets for fresh DNA extracted from HCC1395BL and HCC1395 cell lines",
	"PBO": "PacBio Sequel II Whole Genome Sequencing data sets for fresh DNA extracted from HCC1395BL and HCC1395 cell lines",
	"AMS": "AmpliSeq libraries were prepared using Illumina protocol and sequenced on MiSeq platform",
	"TGEN": "Single cell libraries were prepared using 10X Genomics Chromium Single Cell CNV Solution for CNV profiling",
	"INT": "Ion Torrent"
}

sample_types = {
    "1-0": "100% HCC1395",
    "3-1": "75% HCC1395",
    "1-1": "50% HCC1395",
    "1-4": "25% HCC1395",
    "1-9": "10% HCC1395",
    "1-19": "5% HCC1395",
    "0-1": "0% HCC1395",
    "T": "HCC1395",
    "N": "HCC1395BL"
}

library_prep = {
    "1": "prep #1",
    "2": "prep #2",
    "3": "prep #3",
    "250ng": "250ng of input",
    "100ng": "100ng of input",
    "10ng": "10ng of input",
    "1ng": "1ng of input"
}

sample_annot = (
    meta
    .with_columns(names=pl.col("Library Name").str.split("_"))
    .with_columns(
        study_alias = pl.col("names").list.get(0),
		study_name = pl.col("names").list.get(0).map_elements(lambda x: study_name.get(x, x), return_dtype=str),
		center_alias = pl.col("names").list.get(1),
		sequencing_center = pl.col("names").list.get(1).map_elements(lambda x: sequencing_center.get(x, x), return_dtype=str),
  		sample_types = pl.col("names").list.get(2).map_elements(lambda x: sample_types.get(x, x), return_dtype=str),
    	library_prep = pl.col("names").list.get(3),
	)
    .drop("names")
    .sort("study_alias")
    .select(['study_alias', 'study_name', 'center_alias', 'sequencing_center', 'Instrument', 'sample_types', 'library_prep', 'Tissue',  'Run', 'BioProject', 'BioSample', 'Experiment', 'Library Name', 'Sample Name', 'SRA Study'])
)

sample_annot.write_csv("sample_annotation.tsv", separator="\t")

sample_annot

sample_annot.filter(pl.col("study_alias") == "WGS")


