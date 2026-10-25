cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[2.2745541594922543, 1.1414938680827618, 55.46752538532019, -0.2688746079802513, 2.304817132651806]' OUTTAG=s2_exp_back_gtc_4 timeout 2400 python3 cmb4.py exp back gtc 30 > gr2/v4/s2/s2_exp_back_gtc_4.log 2>&1
