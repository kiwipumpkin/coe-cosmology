cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[50.79738095106633, -1.937881243797385, 2.419649514780748, 3.096887050270368, 0.14863012225111216, 0.6908738454770247]' OUTTAG=s2_pl_fwd_le_3 timeout 2400 python3 cmb4.py pl fwd le 30 > gr2/v4/s2/s2_pl_fwd_le_3.log 2>&1
