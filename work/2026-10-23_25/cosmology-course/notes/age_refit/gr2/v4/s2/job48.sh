cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[2.0455272539021906e-10, 1.1827291565684634, 25.074934246259325, 3.787070874214261, 3.2712491174668434, 5.1377577748301535]' OUTTAG=s2_pl_back_gtc_0 timeout 2400 python3 cmb4.py pl back gtc 30 > gr2/v4/s2/s2_pl_back_gtc_0.log 2>&1
