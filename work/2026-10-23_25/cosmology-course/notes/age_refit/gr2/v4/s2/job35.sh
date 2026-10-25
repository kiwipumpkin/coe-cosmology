cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.6725301565947011, 0.80023183953017, 47.690926399081945, -1.686251476407051, 4.800443686544895]' OUTTAG=s2_exp_back_le_5 timeout 2400 python3 cmb4.py exp back le 30 > gr2/v4/s2/s2_exp_back_le_5.log 2>&1
