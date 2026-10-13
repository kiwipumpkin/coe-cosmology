#!/bin/bash
# relaunch a phase until it reports done (other jobs on the box sometimes kill python processes)
cd /workspace/clump_scan/joint_network
for ph in "$@"; do
  for try in $(seq 1 30); do
    python3 run_jn.py 7 $ph >> log_$ph.txt 2>&1 < /dev/null
    grep -q "^done" log_$ph.txt && break
    sleep 5
  done
done
echo ALLDONE >> log_keep.txt
