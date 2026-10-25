cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[1.0780077809731827e-10, 1.1827282614463732, 44.00725050406992, -2.247491180030039, 4.38069603266261]' OUTTAG=s2_exp_back_gtc_0 timeout 2400 python3 cmb4.py exp back gtc 30 > gr2/v4/s2/s2_exp_back_gtc_0.log 2>&1
