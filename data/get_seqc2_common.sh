#!/bin/bash

# View files to download
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG

# Download files
wget -c -r -nH --cut-dirs=3 --no-parent --execute robots=off "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/assembly"

wget -c -r -nH --cut-dirs=4 --no-parent --execute robots=off "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/technical"

wget -c -r -nH --cut-dirs=3 --no-parent --execute robots=off "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release"

wget -c -r -nH --cut-dirs=3 --no-parent --execute robots=off "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/tools"

mv technical/reference_genome ./seqc2-reference-genome
