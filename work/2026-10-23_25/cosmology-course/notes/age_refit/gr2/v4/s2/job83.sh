cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.9348030093745177, 0.885802299939827, 39.40263836412373, 1.261422795585099, 3.55248938657284, 3.4904511950566453]' OUTTAG=s2_soft_back_le_5 timeout 2400 python3 cmb4.py soft back le 30 > gr2/v4/s2/s2_soft_back_le_5.log 2>&1
