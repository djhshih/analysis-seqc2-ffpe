#!/bin/bash

## This script filters the mutect2 variants after the FilterMutectCalls stage. 
## Variants with "PASS" and "Orientation" annotation in the Filter Column are retained

set -euo pipefail

filter_expression="FILTER='PASS' | FILTER='orientation'"

for vcf in ../mutect2-matched-normal_filtermutectcalls_obmm_unfiltered/*/*/*.vcf; do

    echo -e "\nProcessing $vcf \n"

    base=$(basename $vcf)
    dir_name=$(dirname $vcf)
    file_name=${base%.*}

    sample_name=$(basename $dir_name)
    study_name=$(basename $(dirname $dir_name))

    out_dir="${study_name}/${sample_name}"
    mkdir -p $out_dir

    out_path=${out_dir}/${file_name}.vcf

    bcftools view -i "$filter_expression" $vcf -o $out_path

    echo -e "\nFiltered vcf saved at: $out_path \n"

    gatk IndexFeatureFile -I $out_path

    echo -e "\nFinished Indexing $out_path \n"

done
