cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.83931958823206, -1.6440457600718612, 1.023686931626232, 0.9853725902637762, -2.363464407143219, 3.555038521611453]' OUTTAG=s2_soft_fwd_le_1 timeout 2400 python3 cmb4.py soft fwd le 30 > gr2/v4/s2/s2_soft_fwd_le_1.log 2>&1
