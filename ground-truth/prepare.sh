#!/usr/bin/env bash

set -euo pipefail


cd individual
python create-ground-truth.py

cd ../wes-wgs-union
python create-ground-truth.py

cd ../from-seqc
python create-ground-truth.py

