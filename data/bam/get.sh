#!/bin/bash

set -euo pipefail

# View files to download
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG

# Download files
# wget -c -r -nH --cut-dirs=4 --no-parent --execute robots=off "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data"

bash get_wgs.sh
bash get_wes.sh
bash get_ffg.sh
bash get_ffx.sh
bash get_lbp.sh
bash get_ams.sh
