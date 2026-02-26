#!/usr/bin/env bash

set -euox pipefail

bash filter_pass-orientation.sh
bash wes_filter_exome.sh
bash wes_filter_dp20.sh
bash wes_filter_blacklist.sh
bash wes_filter_clonal.sh
bash wgs_filter_dp20.sh
bash wgs_filter_blacklist.sh
bash wgs_filter_clonal.sh

