#!/bin/bash

## This script filters the variants only retaining those in the exome target region

set -euo pipefail

regions="../../data/blacklists/master_blacklist.bed"

for vcf in ../mutect2-matched-normal_hc-target_pass-orientation-dp-filtered/{WGS,FFG}/*/*.vcf; do

    echo -e "\nProcessing $vcf \n"

    base=$(basename $vcf)
    dir_name=$(dirname $vcf)
    file_name=${base%.*}

    sample_name=$(basename $dir_name)
    study_name=$(basename $(dirname $dir_name))

    out_dir="${study_name}/${sample_name}"
    mkdir -p $out_dir

    out_path=${out_dir}/${file_name}.vcf

    if [[ -f $out_path ]]; then
        echo $out_path already exists.
        continue
    fi

    bcftools view -T ^"$regions" $vcf -o $out_path

    echo -e "\nFiltered vcf saved at: $out_path \n"

    gatk IndexFeatureFile -I $out_path

    echo -e "\nFinished Indexing $out_path \n"

done
