# Analysis SEQC2 FFPE

Here we use the datasets with matched FFPE and Fresh Frozen (FF) samples from the SEQC2 study to benchmark the performance of FFPE artifact filtering tools against our in house orientation-bias based statistical model - MOBSNVF.

## Study

[Whole genome and exome sequencing reference datasets from a multi-center and cross-platform benchmark study](https://www.nature.com/articles/s41597-021-01077-5)

The study generated a comprehensive set of WGS and WES data using two well-characterized reference samples (paired tumor-normal): a human breast cancer cell line (HCC1395) and a B lymphocytes cell line (HCC1395BL) derived from the same donor.

WGS and WES data were generated using various NGS library preparation protocols, seven NGS platforms (NovaSeq, HiSeq, PacBio, 10X Genomics, Ion Torrent, Miseq, and Affymetrix CytoScan HD) at six centers including Illumina (IL), National Cancer Institute (NC), Novartis (NV), European Infrastructure for Translational Medicine (EA), Fudan University (FD), and Loma Linda University (LL).

DNA was extracted from fresh cells or cell pellets mimicking the formalin-fixed paraffin-embedded (FFPE) process with fixation time of 1, 2, 6, or 24 hours. A small amount of DNA from fresh cells of HCC1395 and HCC1395BL was pooled at various ratios (3:1, 1:1, 1:4, 1:9 and 1:19) to create mixtures. Both fresh DNA and FFPE DNA were profiled on NGS or microarray platforms following manufacturer recommended protocols. To assess the reproducibility of WGS and WES, six sequencing centers performed a total of 12 replicates (3 × 3 + 3) on each platform. In addition, 12 WGS libraries constructed using three different library preparation protocols (TruSeq PCR-free, TruSeq-Nano, and Nextera Flex) in four different quantities of DNA inputs (1, 10, 100, and 250 ng) were sequenced on an Illumina HiSeq 4000, and nine WGS libraries constructed using the TruSeq PCR-free protocol were sequenced on an Illumina NovaSeq. Finally, Affymetrix Cytoscan HD and single-cell sequencing with 10X Genomics platform were performed to uncover the cytogenetics and heterogeneity of two cell lines.

![Zhao et. al. 2021](image.png)

## Datasets

The study has a dedicated webpage at: https://sites.google.com/view/seqc2/home/
Additionally, the data is deposited at: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/

- **WGS:** Whole Genome Sequencing (WGS) data sets for fresh DNA extracted from HCC1395BL and HCC1395 cell lines

- **WES:** Whole Exome Sequencing (WES) data sets for fresh DNA extracted from HCC1395BL and HCC1395 cell lines

- **FFG:** Whole Genome Sequencing (WGS) data sets for Formalin-Fixed Paraffin-Embedded (FFPE) process with fixation time of 1, 2, 6, or 24 hours for DNA extracted from HCC1395BL and HCC1395 cell lines

- **FFX:** Whole Exome Sequencing (WES) data sets for DNA extracted from HCC1395BL and HCC1395 cell lines and processed via Formalin-Fixed Paraffin-Embedded (FFPE) process with fixation time of 1, 2, 6, or 24 hours

- **LBP:** WGS Libraries were made from different library protocols such as TruSeq Nano, TruSeq PCR Free and Nextera Flex library protocol with different input amount and sequenced on Illumina HiSeq 3000/4000

- **SPP:** WGS Libraries were made from pooling the HCC1395 and HCC1395BL cell lines with various ratios (3:1, 1:1, 1:4, 1:9 and 1:19) to create mixtures

- **AMS:** AmpliSeq libraries were prepared using Illumina protocol and sequenced on MiSeq platform.

## Accession

- **SRA:** [SRP162370](https://www.ncbi.nlm.nih.gov/sra/?term=SRP162370)
- **ENA:** [SRP162370](https://www.ebi.ac.uk/ena/browser/view/SRP162370)
- **NCBI FTP:** https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG


## Cell lines

Cell lines HCC1395 (tumour) and HCC1395BL (blood) were derived from the same patient
with stage 1 breast ductal carcinoma.

HCC1395
https://www.atcc.org/products/crl-2324

HCC1395BL
https://www.atcc.org/products/crl-2325
