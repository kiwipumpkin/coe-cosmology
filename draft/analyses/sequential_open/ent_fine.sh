#!/bin/bash
# stage 2 (fine-scale pass) for every web whose stage-1 run is finished; resumes from ckpt
cd /workspace/clump_scan/sequential_open
for t in "$@"; do
  mv -f minedit_ent/${t}_tidalcore_entropy.json minedit_ent/stage1_${t}.json 2>/dev/null
  FINE=1 OMP_NUM_THREADS=1 python3 minedit_ent.py ${t}_tidalcore > minedit_ent/log2_${t}.txt 2>&1 &
done
wait
