# Phase5 surface-definition test (Oct 26 PM): background-only, SURF in {t,rho} for K-1, A_d2, poly validation. nohup; resume-safe (bg_mig skips done traces).
cd /workspace/coepdf/joint_solve/phase5
run(){ T=$1;S=$2;shift 2; [ -f bg_${T}_ckpt.json ] || true; [ -f bg_$T.json ] || SURF=$S OMP_NUM_THREADS=1 python3 bg_mig.py $T "$@" > surf/$T.log 2>&1; }
for S in t rho; do
 cp -n bg_K-1_ckpt.json bg_K-1_S$S\_ckpt.json; cp -n bg_A_d2_ckpt.json bg_A_d2_S$S\_ckpt.json
 run K-1_S$S $S 27.737393 0.060341 3.068065 11.731858 &
 run A_d2_S$S $S 27.354588 0.056062 9.127586 12.425732 &
done
run VALpolyA_St t 27.3385 0.0185 1 1 poly &
wait; echo SURF_ALL_DONE
