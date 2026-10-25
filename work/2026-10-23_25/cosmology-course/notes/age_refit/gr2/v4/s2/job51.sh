cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.6507397007380025, 0.8770570763124186, 26.572478080409162, 4.67532741101095, 4.331211307223177, 5.16192034857004]' OUTTAG=s2_pl_back_gtc_3 timeout 2400 python3 cmb4.py pl back gtc 30 > gr2/v4/s2/s2_pl_back_gtc_3.log 2>&1
