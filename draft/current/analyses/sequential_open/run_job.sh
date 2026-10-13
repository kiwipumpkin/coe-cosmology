#!/bin/bash
# usage: run_job.sh MODE KIND RULE SEED   (skips if output exists)
cd /workspace/clump_scan/sequential_open
m=$1; k=$2; r=$3; s=$4; suf=""; [ "$m" != both ] && suf="_$m"
export SEQ_OPEN=$m OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
if [ "$k" = b ]; then out=b_${r}_s${s}${suf}.json; [ -f $out ] || python3 -u run_b_open.py $r $s >> log_b_${r}_s${s}${suf/_/_}.txt 2>&1
else out=shallow_${r}_s${s}${suf}.json; [ -f $out ] || python3 -u run_shallow_open.py grow $r $s >> log_shallow_${r}_s${s}${suf}.txt 2>&1; fi
