cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[50.60487808642985, -2.4705917224384777, 1.7542603318033696, 3.141760669722487, 1.302564834772844, 0.8901192648302714]' OUTTAG=s2_pl_fwd_le_0 timeout 2400 python3 cmb4.py pl fwd le 30 > gr2/v4/s2/s2_pl_fwd_le_0.log 2>&1
