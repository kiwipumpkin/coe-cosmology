cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[54.129447623662884, -2.0803396724353465, 1.2081129024637383, 2.24195753721059, 5.973988802180459, 5.2610787290665435]' OUTTAG=s2_soft_fwd_gtc_4 timeout 2400 python3 cmb4.py soft fwd gtc 30 > gr2/v4/s2/s2_soft_fwd_gtc_4.log 2>&1
