#!/bin/bash
# bodies-on-top with added bodies kept >= RING_EX Gpc (surface) from the CMB ring; usage: omp_queue_ex.sh TAG...
cd /workspace/clump_scan/sequential_open
for t in "$@"; do [ -f final_$t.npz ] && [ ! -f web_add/${t}_omp0_ex1.5.json ] && echo $t; done | xargs -P ${P:-4} -I{} sh -c 'RING_EX=1.5 OSUF=_ex1.5 NBANK=80000 OMP_NUM_THREADS=2 python3 -u omp_add.py {} 0 >> log_omp_ex.txt 2>&1'
