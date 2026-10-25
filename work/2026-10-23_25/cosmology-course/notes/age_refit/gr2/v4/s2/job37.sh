cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[55.68541427092013, -1.8976169172004438, 1.4551570204979458, -0.7359881883090065, 0.34017547321533614]' OUTTAG=s2_exp_fwd_gtc_1 timeout 2400 python3 cmb4.py exp fwd gtc 30 > gr2/v4/s2/s2_exp_fwd_gtc_1.log 2>&1
