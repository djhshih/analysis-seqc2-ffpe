#!/usr/bin/env bash

wget -c -r -nH --cut-dirs=6 --no-parent --execute robots=off "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/"

