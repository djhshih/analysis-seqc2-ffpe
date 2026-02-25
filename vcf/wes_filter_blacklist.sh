#!/usr/bin/env bash

set -euo pipefail

blacklist_regions="../data/blacklists/master_blacklist.bed"
filtered_set_name="mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist"

i=1

for vcf in {WES,FFX}/mutect2-tn_filtered_pass-orientation-exome-dp20/*/*.vcf.gz; do

    echo "$i. Filtering: $vcf"

    dataset=$(echo $vcf | cut -d'/' -f1)
    sample_name=$(echo $vcf | cut -d'/' -f3)

    outdir=$dataset/$filtered_set_name/$sample_name
    mkdir -p $outdir
    
    outpath=${outdir}/${sample_name}.vcf.gz

    bcftools view -T "$blacklist_regions" -Oz -o $outpath $vcf

    echo -e "\t Indexing..."
    bcftools index -t $outpath

    echo -e "\t Filtered vcf saved at: $outpath \n"

    ((i++))

done
