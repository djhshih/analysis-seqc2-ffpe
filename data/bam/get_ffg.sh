#!/bin/bash

set -euo pipefail

# Download files
mkdir FFG
cd FFG
wget -c -r -nH --cut-dirs=4 --no-parent --execute robots=off "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/FFG"
mv FFG no_rg_bam
