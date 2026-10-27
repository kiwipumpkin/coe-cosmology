# Phase5 routes (NEW): background-only builds (bg_mig.py) for each route's best fit, SURF=t (primary) and SURF=a. Resume-safe. Usage: bash routes/run_bg.sh TAG...
cd /workspace/coepdf/joint_solve/phase5
for T in "$@"; do P=$(python3 -c "import json;r=json.load(open('routes/fit_$T.json'));print(r['D'],r['mc'],r['eta'],r['L'])")
 for S in t a; do [ -f bg_R${T}_S$S.json ] || SURF=$S OMP_NUM_THREADS=1 python3 bg_mig.py R${T}_S$S $P > routes/bg_R${T}_S$S.log 2>&1 & done; done; wait; echo BG_DONE
