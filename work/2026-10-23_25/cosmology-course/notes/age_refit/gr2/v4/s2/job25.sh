cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.31497730325511686, 1.2003368742221758, 44.47550087837964, -1.0462937346374304, 4.77736166818134]' OUTTAG=s2_exp_back_gtc_1 timeout 2400 python3 cmb4.py exp back gtc 30 > gr2/v4/s2/s2_exp_back_gtc_1.log 2>&1
