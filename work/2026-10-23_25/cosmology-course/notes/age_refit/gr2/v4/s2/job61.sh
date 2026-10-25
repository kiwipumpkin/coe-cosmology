cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[54.535017326987145, -1.987776212975049, 0.46901065484352755, 0.3838736279978532, 0.44462717616480907, 3.7149304017195046]' OUTTAG=s2_pl_fwd_gtc_1 timeout 2400 python3 cmb4.py pl fwd gtc 30 > gr2/v4/s2/s2_pl_fwd_gtc_1.log 2>&1
