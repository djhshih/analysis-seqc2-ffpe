#!/usr/bin/env bash

set -euo pipefail

filter_expression='MIN(FMT/DP)>=20'
source_set_name="mutect2-tn_filtered_pass-orientation"
filtered_set_name="mutect2-tn_filtered_pass-orientation-dp20"

i=1

for vcf in {WGS,FFG}/$source_set_name/*/*.vcf.gz; do

    echo "$i. Filtering: $vcf"

    dataset=$(echo $vcf | cut -d'/' -f1)
    sample_name=$(echo $vcf | cut -d'/' -f3)

    outdir=$dataset/$filtered_set_name/$sample_name
    mkdir -p $outdir
    
    outpath=${outdir}/${sample_name}.vcf.gz

    bcftools view -i "$filter_expression" -Oz -o $outpath $vcf

    echo -e "\t Indexing..."
    bcftools index -t $outpath

    echo -e "\t Filtered vcf saved at: $outpath \n"

    ((i++))

done
