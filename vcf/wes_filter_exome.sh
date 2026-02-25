#!/usr/bin/env bash

set -euo pipefail

exome_regions="../data/seqc2-reference-genome/Exome_Target_bed/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed"
filtered_set_name="mutect2-tn_filtered_pass-orientation-exome"

i=1

for vcf in {WES,FFX}/mutect2-tn_filtered_pass-orientation/*/*.vcf.gz; do

    echo "$i. Filtering: $vcf"

    dataset=$(echo $vcf | cut -d'/' -f1)
    sample_name=$(echo $vcf | cut -d'/' -f3)

    outdir=$dataset/$filtered_set_name/$sample_name
    mkdir -p $outdir
    
    outpath=${outdir}/${sample_name}.vcf.gz

    bcftools view -T "$exome_regions" -Oz -o $outpath $vcf

    echo -e "\t Indexing..."
    bcftools index -t $outpath

    echo -e "\t Filtered vcf saved at: $outpath \n"

    ((i++))

done
