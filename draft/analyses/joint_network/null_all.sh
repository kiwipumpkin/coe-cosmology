#!/bin/bash
cd /workspace/clump_scan/joint_network
for c in "cont m 12" "cont n 20" "web m 24" "loose m 24" "strict m 24" "web n 24"; do
  set -- $c; f=null/${1}_${2}_N${3}.json
  [ -f $f ] && continue
  for try in 1 2 3; do python3 null_nets.py 7 100 $1 $2 $3 >> log_null.txt 2>&1 < /dev/null && break; done
done
echo NULLDONE >> log_null.txt
