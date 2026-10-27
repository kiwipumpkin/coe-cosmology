# Sequential (OOM history) inverse runs; resume-safe via per-seed checkpoints; MAXMV=1000 is a COMPUTE budget (1/3 of bodies), not physics.
cd /workspace/coepdf/joint_solve/phase5/web/inverse; export OMP_NUM_THREADS=4
for s in 1 2 3; do [ -f inv_s$s.json ] || python3 inverse2.py --seed $s --MAXMV 1000 >> inv_s$s.log 2>&1 || { echo FAIL $s; exit 1; }; done
python3 tt_inverse.py 1 2 3 > tt_inverse.log 2>&1 && python3 TT_inverse_log.py && python3 radial_moved.py && echo ALL_DONE
