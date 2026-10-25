#!/bin/bash
# runO.sh A LAM LD PHTAG(trough|crest|up|down) [frac]
f=${5:-30.0}; n=$PWD/runs_osc/o_$1_$2_$3_$4_f$f.npz; [ -f $n ] && exit 0
PH=$(python3 -c "import math;l=$2;d={'trough':-math.pi/2,'crest':math.pi/2,'up':0,'down':math.pi}['$4'];print(d-2*math.pi/l)")
A=$1 LAM=$2 LD=$3 PH=$PH python3 famO.py osc 0 $f $n >> osc.log 2>&1
