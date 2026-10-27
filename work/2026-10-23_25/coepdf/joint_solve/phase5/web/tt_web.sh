# Phase5 web test: full TT chain on Route A eps=0.01 St background after web compensation (cap 1e22, 0% of bodies at cap). Resume-safe, fail-loud.
set -x;set -o pipefail;cd /workspace/coepdf/joint_solve/phase5;T=web_RA0.01_M1e22;fail(){ echo ONE_FAIL $T "$1"; exit 1; }
D=30.424559192759354;R=21.7439;C="--bg $PWD/bg_RA_0.01_St_sm_ns64.npy --R $R --MMAX 1e22 --LC 40 --LP 60 --WH 1 --NB 3000 --obs 1 --LQ 20 --LT 40"
python3 check_finite.py p5_$T.npz || fail p5_nonfinite
[ -f p8_$T.json ] || TAG=$T DG=$D R=$R DFIXF=$PWD/p5_$T.npz NP=${NP:-6} python3 p8_scratch5.py || fail p8
python3 check_finite.py p8_$T.npz || fail p8_nonfinite
[ -f p5_fld_$T.json ] || python3 p5_comp5.py $C --field $PWD/p8_$T.npz --f 1 --tag fld_$T || fail p5_fld
python3 check_finite.py p5_fld_$T.npz || fail p5_fld_nonfinite
[ -f p6_$T.json ] || python3 p6_tt5.py $PWD/p8_$T.npz $PWD/p5_fld_$T.npz $T || fail p6
python3 check_finite.py p6_$T.npz || fail p6_nonfinite
python3 web/tt_web_summary.py && python3 web/TT_web_plot.py && LOG=1 python3 web/TT_web_plot_log.py && echo ONE_DONE $T
