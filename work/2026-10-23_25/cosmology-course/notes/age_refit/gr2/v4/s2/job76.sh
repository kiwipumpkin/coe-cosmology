cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.17789709238648999, 1.1209323002113172, 46.433160444299574, -2.812295439128472, -2.795878183461549, 5.749871654600886]' OUTTAG=s2_soft_back_gtc_4 timeout 2400 python3 cmb4.py soft back gtc 30 > gr2/v4/s2/s2_soft_back_gtc_4.log 2>&1
