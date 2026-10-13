#!/bin/bash
# detached; each run resumes from ckpt/minedit_ent_<tag>_entropy.pkl if interrupted
cd /workspace/clump_scan/sequential_open
for t in b_random_s1 b_random_s2 b_random_s3 b_foam_rim_s1 b_foam_rim_s2 b_foam_rim_s3; do
  [ -f minedit_ent/${t}_tidalcore_entropy.json ] && continue
  OMP_NUM_THREADS=1 python3 minedit_ent.py ${t}_tidalcore > minedit_ent/log_${t}.txt 2>&1 &
done
wait
