#!/usr/bin/env bash

set -euo pipefail

## This script fixes the amplicon file names to align with the annotation table and adds RG to the header of the FFG bams

cd AMS
python fix_name.py

cd -

cd FFG
bash add_rg.sh