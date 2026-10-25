cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/best XOVR='[0.6507397007380025, 0.8770570763124186, 26.572478080409162, 4.67532741101095, 4.331211307223177, 5.16192034857004]' OUTTAG=best_pl_back_gtc_3 python3 cmb4.py pl back gtc 6 > gr2/v4/best/best_pl_back_gtc_3.log 2>&1
