#!/bin/bash

mkdir -p gatk-test-data/mutect2/
cd gatk-test-data/mutect2/
gsutil cp -nr 'gs://gatk-test-data/mutect2/*' .

cd -
mkdir -p gatk-best-practices/somatic-hg38/
cd gatk-best-practices/somatic-hg38/
gsutil cp -n 'gs://gatk-best-practices/somatic-hg38/*' .

cd -
mkdir -p gatk-reference-genome
cd gatk-reference-genome
gsutil cp -n 'gs://genomics-public-data/resources/broad/hg38/v0/*' .