# Seed robustness: same as run_web.sh (eps=0.01, cap 1e22) but --seed varied; sequential (memory); skips finished seeds (checkpoint).
cd /workspace/coepdf/joint_solve/phase5; export OMP_NUM_THREADS=2
R=$(python3 -c "import json;print(round(json.load(open('bg_RA_0.01_St.json'))['R_ls'],4))")
for S in 11 23 101 777 1234 3141 5000 9001; do T=web_RA0.01_M1e22_s$S; [ -f p5_$T.json ] && continue
 python3 p5_comp5.py --bg $PWD/bg_RA_0.01_St_sm_ns64.npy --R $R --MMAX 1e22 --LC 40 --LP 60 --WH 1 --NB 3000 --obs 1 --LQ 20 --LT 40 --seed $S --tag $T > web/$T.log 2>&1
done; echo SEEDS_DONE
