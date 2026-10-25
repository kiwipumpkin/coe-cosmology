cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.625335337935466, 1.2326976109675922, 38.37452242198083, 1.6767167275461585, 3.946823819241824, 4.427503571624962]' OUTTAG=s2_soft_back_le_1 timeout 2400 python3 cmb4.py soft back le 30 > gr2/v4/s2/s2_soft_back_le_1.log 2>&1
