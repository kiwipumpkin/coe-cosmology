cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[50.037677140739575, -2.318615443455421, 1.4218794819353104, 4.2426031363003744, 2.6648523002557942, 0.30521263833753265]' OUTTAG=s2_pl_fwd_le_2 timeout 2400 python3 cmb4.py pl fwd le 30 > gr2/v4/s2/s2_pl_fwd_le_2.log 2>&1
