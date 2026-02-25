#!/usr/bin/env bash

set -euo pipefail

source_set_name="mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist"
filtered_set_name="mutect2-tn_filtered_pass-orientation-exome-dp20-blacklist-clonal"

clonal_variants_root="../stratify-clonal-subclonal"

i=1

for vcf in {WES,FFX}/$source_set_name/*/*.vcf.gz; do

    echo "$i. Filtering: $vcf"

    dataset=$(echo $vcf | cut -d'/' -f1)
    sample_name=$(echo $vcf | cut -d'/' -f3)

    outdir=$dataset/$filtered_set_name/$sample_name
    mkdir -p $outdir
    
    outpath=${outdir}/${sample_name}.vcf.gz

    clonal_variants=$clonal_variants_root/$dataset/$source_set_name/$sample_name/$sample_name.clonal.tsv

    if [[ ! -f $clonal_variants ]]; then
        echo Clonal variant table not found at: $clonal_variants
        echo execution halted
        exit 1
    fi

    bcftools view -T "$clonal_variants" -Oz -o $outpath $vcf

    echo -e "\t Indexing..."
    bcftools index -t $outpath

    echo -e "\t Filtered vcf saved at: $outpath \n"

    ((i++))

done
