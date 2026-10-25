cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.6521917452251076, 1.140543235225396, 42.80684283791572, -1.1075225901504666, 5.589473009575824]' OUTTAG=s2_exp_back_gtc_2 timeout 2400 python3 cmb4.py exp back gtc 30 > gr2/v4/s2/s2_exp_back_gtc_2.log 2>&1
