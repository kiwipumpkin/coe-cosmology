#!/bin/bash
cd /workspace/clump_scan/solution_space
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
xargs -P 8 -L 1 sh -c 'python3 fit_run.py $0 $1 $2 2>&1 | grep -v beta=0' < jobs2.txt >> log_fits2.txt
echo ALLDONE >> log_fits2.txt
