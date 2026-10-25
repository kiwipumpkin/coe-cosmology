#!/bin/bash
n=$PWD/runs_poly/f$5_p_$1_$2_$3_$4.npz; [ -f $n ] && exit 0
C1=$1 C2=$2 C3=$3 C4=$4 python3 famP.py poly 0 $5 $n >> poly.log 2>&1
