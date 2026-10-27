# Phase5 (NEW): full CMB TT pipeline on the background built from the Migkas density profile (bg_mig.py). Resume-safe (each step skipped if output exists).
# Usage: bash cmb_mig.sh TAG D mc eta L_Gpc
set -x;set -o pipefail;fail(){ echo ONE_FAIL $T "$1"; exit 1; };cd /workspace/coepdf/joint_solve/phase5;T=$1;D=$2;MC=$3;ETA=$4;L=$5;NP=${NP:-3}
[ -f bg_$T.json ] || python3 bg_mig.py $T $D $MC $ETA $L || fail bg_mig
python3 check_finite.py bg_${T}_sm_ns64.npy bg_${T}_trace.npz || fail bg_nonfinite
python3 -c "import json,sys;b=json.load(open('bg_$T.json'));sys.exit(1 if 'fail' in b else 0)" || { echo ONE_FAIL $T; exit 1; }
R=$(python3 -c "import json;print(round(json.load(open('bg_$T.json'))['R_ls'],4))")
C="--bg $PWD/bg_${T}_sm_ns64.npy --R $R --MMAX 1e22 --LC 40 --LP 60 --WH 1 --NB 3000 --obs 1 --LQ 20 --LT 40"
[ -f p5_bg_$T.json ] || python3 p5_comp5.py $C --tag bg_$T || fail p5_bg
python3 check_finite.py p5_bg_$T.npz || fail p5_bg_nonfinite
[ -f p8_$T.json ] || TAG=$T DG=$D R=$R DFIXF=$PWD/p5_bg_$T.npz NP=$NP python3 p8_scratch5.py || fail p8
python3 check_finite.py p8_$T.npz || fail p8_nonfinite
[ -f p5_fld_$T.json ] || python3 p5_comp5.py $C --field $PWD/p8_$T.npz --f 1 --tag fld_$T || fail p5_fld
python3 check_finite.py p5_fld_$T.npz || fail p5_fld_nonfinite
[ -f p6_$T.json ] || python3 p6_tt5.py $PWD/p8_$T.npz $PWD/p5_fld_$T.npz $T || fail p6
python3 check_finite.py p6_$T.npz || fail p6_nonfinite
echo ONE_DONE $T
