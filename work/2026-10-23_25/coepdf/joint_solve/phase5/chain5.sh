# SUPERSEDED 2026-10-26 01:4x CDT: old CMB-from-D-only rebuild disabled (Chris: CMB must use the Migkas density profile). Use cmb_mig.sh / chain6.sh.
echo "chain5 superseded"; exit 0
# Phase5 (NEW) master chain: wait for Migkas fixed-d fits, then CMB pipeline per unique D (rounded 0.01 Gpc), then summary.
cd /workspace/coepdf/joint_solve/phase5
while [ $(ls fitd_*[0-9].json 2>/dev/null | wc -l) -lt 9 ]; do sleep 60; done
python3 -c "
import json,glob
Ds=sorted({round(json.load(open(f))['D']*20)/20 for f in glob.glob('fitd_*[0-9].json')});open('Dlist.txt','w').write('\n'.join('%.2f'%x for x in Ds))"
for D in $(cat Dlist.txt); do bash cmb_one.sh $D D$D; done
python3 summary5.py; python3 chi2_vs_d_plot.py; python3 TT_best_plot.py; python3 TT_best_plot_log.py
echo CHAIN5_DONE
