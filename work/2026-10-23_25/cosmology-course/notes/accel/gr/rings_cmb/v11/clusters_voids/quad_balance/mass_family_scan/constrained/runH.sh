#!/bin/bash
n=$PWD/runs_poly/hyb_xs$1_f$2.npz; [ -f $n ] && exit 0
XS=$1 python3 famP.py hyb 0 $2 $n >> hyb.log 2>&1
