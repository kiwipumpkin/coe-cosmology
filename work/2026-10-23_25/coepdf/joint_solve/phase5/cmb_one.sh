# Phase5 (NEW): CMB pipeline for one background D (Gpc). Each step skipped if its output exists (resume-safe).
# Usage: bash cmb_one.sh D TAG
set -x;cd /workspace/coepdf/joint_solve/phase5;D=$1;T=$2;FR=0.0185
U=$(python3 -c "print(68.64*$D*1000/299792.458)")
[ -f ltb_scan_$T.json ] || DSHAPE=$D TAG=$T python3 ltb5.py "$U,$FR"
R=$(python3 -c "import json;print(round(json.load(open('ltb_scan_$T.json'))[0]['eta_rec_D']*$D,4))")
BG=$(python3 -c "print('p2_sm_D%g_H68.64_f0.0185_ns64.npy'%$D)")
[ -f $BG ] || python3 p2_bg5.py $D 68.64 $FR
C="--bg $PWD/$BG --R $R --MMAX 1e22 --LC 40 --LP 60 --WH 1 --NB 3000 --obs 1 --LQ 20 --LT 40"
[ -f p5_bg_$T.json ] || python3 p5_comp5.py $C --tag bg_$T
[ -f p8_$T.json ] || TAG=$T DG=$D R=$R DFIXF=$PWD/p5_bg_$T.npz NP=6 python3 p8_scratch5.py
[ -f p5_fld_$T.json ] || python3 p5_comp5.py $C --field $PWD/p8_$T.npz --f 1 --tag fld_$T
[ -f p6_$T.json ] || python3 p6_tt5.py $PWD/p8_$T.npz $PWD/p5_fld_$T.npz $T
echo ONE_DONE $T
