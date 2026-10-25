cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[9.422024600207806, 1.212935235351324, 51.62184553220868, -0.3913499116897583, 1.7708861157298088, 2.9558838246390224]' OUTTAG=s2_pl_back_gtc_1 timeout 2400 python3 cmb4.py pl back gtc 30 > gr2/v4/s2/s2_pl_back_gtc_1.log 2>&1
