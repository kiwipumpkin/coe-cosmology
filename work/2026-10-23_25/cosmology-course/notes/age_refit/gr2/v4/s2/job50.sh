cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.34127202929895084, 1.2239523873401839, 25.54436213826626, 5.0906213429720095, 4.725545739892161, 6.098972725138357]' OUTTAG=s2_pl_back_gtc_2 timeout 2400 python3 cmb4.py pl back gtc 30 > gr2/v4/s2/s2_pl_back_gtc_2.log 2>&1
