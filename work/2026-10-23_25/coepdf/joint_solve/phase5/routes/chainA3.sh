# waits for fit_A_0.001.json then builds SURF=t/a backgrounds and metrics (resume-safe)
cd /workspace/coepdf/joint_solve/phase5
while [ ! -f routes/fit_A_0.001.json ]; do sleep 60; done
bash routes/run_bg.sh A_0.001
python3 surf/metrics.py RA_0.001_St RA_0.001_Sa > routes/metrics_RA_0.001.json 2>&1; echo CHAIN_A3_DONE
