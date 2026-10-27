# Sequential v3 chain (OOM history): waits for v2 run_inv.sh (pid 1950803) and inverse2_replay.py (pid 1986187) to exit, then one job at a time.
# Resume-safe: each seed skips if inv3_s<seed>.json exists, else resumes from its checkpoint. TT + plots refreshed after each seed.
cd /workspace/coepdf/joint_solve/phase5/web/inverse/v3; export OMP_NUM_THREADS=6
while kill -0 1950803 2>/dev/null || kill -0 1986187 2>/dev/null; do sleep 60; done; echo "v2 jobs done $(date)"
for s in 1 2 3; do
  [ -f inv3_s$s.json ] || python3 inverse3.py --seed $s >> inv3_s$s.log 2>&1 || { echo "FAIL seed $s $(date)"; exit 1; }
  ss=$(ls inv3_s*.json | sed 's/inv3_s\([0-9]*\).json/\1/' | tr '\n' ' ')
  python3 tt_inverse3.py $ss > tt_inverse3.log 2>&1 && python3 rms_vs_moves.py && python3 radial.py && python3 moved_map.py 1 && python3 TT_log.py || { echo "FAIL post $s"; exit 1; }
  echo "seed $s done $(date)"
done; echo ALL_DONE
