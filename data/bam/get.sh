#!/bin/bash

# View files to download
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG

# Download files
wget -c -r -nH --cut-dirs=4 --no-parent --execute robots=off "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data"
