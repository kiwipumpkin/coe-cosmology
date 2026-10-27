cd /workspace/coepdf/joint_solve/phase5
for d in 0.001 0.5 1 2 3 4 6 8 10; do [ -f fitd_$d.json ] || nohup python3 fitd5.py $d > fitd_$d.log 2>&1 & done; wait; echo FITS_DONE
