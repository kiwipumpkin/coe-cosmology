cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.2806815087010621, 0.8757634095305369, 43.6935283601889, -3.0338167278742105, 3.11397883156381]' OUTTAG=s2_exp_back_gtc_3 timeout 2400 python3 cmb4.py exp back gtc 30 > gr2/v4/s2/s2_exp_back_gtc_3.log 2>&1
