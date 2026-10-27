# Phase5 web test: stochastic-web compensation bodies (p5_comp5.py, same settings as cmb_mig.sh) on Route A St backgrounds, cap 1e22 (OUR bound) vs relaxed.
cd /workspace/coepdf/joint_solve/phase5; export OMP_NUM_THREADS=2
for E in 0.01 0.1; do R=$(python3 -c "import json;print(round(json.load(open('bg_RA_${E}_St.json'))['R_ls'],4))")
 for M in 1e22 1e24 1e27; do T=web_RA${E}_M$M; [ -f p5_$T.json ] && continue
  python3 p5_comp5.py --bg $PWD/bg_RA_${E}_St_sm_ns64.npy --R $R --MMAX $M --LC 40 --LP 60 --WH 1 --NB 3000 --obs 1 --LQ 20 --LT 40 --tag $T > web/$T.log 2>&1
 done; done; wait; echo WEB_ALL_DONE
