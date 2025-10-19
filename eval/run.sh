#!/usr/env/bin bash

set -euo pipefail

Rscript eval-gatk-obmm_precrec.R
Rscript eval-mobsnvf_precrec.R
Rscript eval-vafsnvf_precrec.R
Rscript eval-sobdetector_precrec.R

Rscript combine_results.R
Rscript make-plots.R