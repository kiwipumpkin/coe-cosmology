#!/bin/bash
# order: no background and shallow first (all lambdas), then the rest
cd /workspace/clump_scan/sequential_open
cat <<J | xargs -P ${P:-4} -L 1 sh -c 'LAM=$1 OMP_NUM_THREADS=2 python3 -u prior_add.py $0 >> log_prior.txt 2>&1'
none_s0 1
shallow_cmb_aware_s2_tidalcore 1
none_s0 0
shallow_cmb_aware_s2_tidalcore 0
b_random_s2_tidalcore 1
b_foam_rim_s3_tidalcore 1
wave_s1 1
lattice_s1 1
none_s0 20
shallow_cmb_aware_s2_tidalcore 20
b_random_s2_tidalcore 0
b_foam_rim_s3_tidalcore 0
wave_s1 0
lattice_s1 0
b_random_s2_tidalcore 20
b_foam_rim_s3_tidalcore 20
wave_s1 20
lattice_s1 20
J
