#!/bin/bash
# first-pass (drop rule) web-add on the filled webs used in the variant figures + extra seeds; skips finished
cd /workspace/clump_scan/sequential_open
for t in b_random_s2 b_foam_rim_s3 b_random_s2_core b_foam_rim_s3_core b_random_s2_inside b_foam_rim_s3_inside b_cmb_aware_s1 b_cmb_aware_s1_core shallow_cmb_aware_s2 shallow_cmb_aware_s2_core shallow_cmb_aware_s2_inside b_random_s1 b_foam_rim_s1 b_random_s3 b_random_s4 b_random_s5 b_foam_rim_s2; do
  [ -f final_$t.npz ] && [ ! -f web_add/${t}_omp0.json ] && echo $t; done | xargs -P 3 -I{} sh -c 'NBANK=80000 python3 -u omp_add.py {} 0 >> log_omp.txt 2>&1'
