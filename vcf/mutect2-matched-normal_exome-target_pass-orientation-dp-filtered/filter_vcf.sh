#!/bin/bash

## This script filters the variants only retaining those in the exome target region

set -euo pipefail

# filter_expression="(FILTER='PASS' | FILTER='orientation') & MIN(FMT/DP)>=10"

regions="../../data/seqc2-reference-genome/Exome_Target_bed/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed"

for vcf in ../mutect2-matched-normal_pass-orientation-dp-filtered/{WES,FFX}/*/*.vcf; do

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

    bcftools view -T "$regions" $vcf -o $out_path

    echo -e "\nFiltered vcf saved at: $out_path \n"

    gatk IndexFeatureFile -I $out_path

    echo -e "\nFinished Indexing $out_path \n"

done
