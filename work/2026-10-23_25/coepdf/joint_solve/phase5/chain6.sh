# Phase5 (NEW) master chain: CMB on Migkas-profile backgrounds (SURF=a), 3 in parallel (NP=2 each), resume-safe. Re-queues once all fitdk fits exist.
cd /workspace/coepdf/joint_solve/phase5;export SURF=a
run(){ python3 queue6.py > queue6.txt; xargs -P 3 -L 1 bash -c 'NP=2 bash cmb_mig.sh "$0" "$@" > cmb_$0.log 2>&1' < queue6.txt; }
run
while [ $(ls fitdk_*[0-9].json 2>/dev/null | wc -l) -lt 9 ]; do sleep 120; done
run
grep -l ONE_FAIL cmb_*.log && echo CHAIN6_HAS_FAILURES;python3 summary6.py; python3 chi2_vs_d6_plot.py; python3 TT_best6_plot.py; python3 TT_best6_plot_log.py
echo CHAIN6_DONE
