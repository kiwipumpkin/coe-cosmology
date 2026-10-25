cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.17789707725833814, 1.1209349780080873, 25.199377043413513, 3.512291058056145, 3.2704810026541713, 4.638563760012677]' OUTTAG=s2_pl_back_gtc_4 timeout 2400 python3 cmb4.py pl back gtc 30 > gr2/v4/s2/s2_pl_back_gtc_4.log 2>&1
