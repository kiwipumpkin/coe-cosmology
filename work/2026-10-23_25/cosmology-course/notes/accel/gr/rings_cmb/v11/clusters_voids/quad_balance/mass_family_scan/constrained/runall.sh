#!/bin/bash
r(){ if [ "$1" = H ]; then bash runH.sh $2 $3; else bash runP2.sh "$@"; fi; }
export -f r
cat interp_jobs.txt hyb_jobs.txt | xargs -P 7 -L 1 bash -c 'r "$@"' _
echo ALLDONE
