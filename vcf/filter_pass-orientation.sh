#!/usr/bin/env bash

set -euo pipefail

filter_expression="FILTER='PASS' | FILTER='orientation'"
filtered_set_name="mutect2-tn_filtered_pass-orientation"

i=1

for vcf in {AMS,LBP,WES,WGS,FFX,FFG}/mutect2-tumor-normal_filtermutectcalls_obmm_unfiltered/*/*.vcf; do

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
