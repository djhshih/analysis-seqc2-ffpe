#!/bin/bash

set -euo pipefail

# Download files
wget -c -r -nH --cut-dirs=4 --no-parent --execute robots=off "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS"