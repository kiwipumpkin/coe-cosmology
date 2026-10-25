cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.46209540293738244, 1.2225275672972202, 46.22650919482112, 3.0757324025034904, 0.6272245943546295]' OUTTAG=s2_exp_back_le_1 timeout 2400 python3 cmb4.py exp back le 30 > gr2/v4/s2/s2_exp_back_le_1.log 2>&1
