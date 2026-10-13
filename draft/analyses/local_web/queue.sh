#!/bin/bash
# sequential job queue: queue.sh "<args1>" "<args2>" ...   (each args = s nreal seed mode tag)
cd /workspace/clump_scan/local_web
for a in "$@"; do set -- $a; python3 run_web.py $1 $2 $3 $4 $5 $6 > log_$5.txt 2>&1; done
