cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.2840633296538337, 1.1914734450372575, 37.934816416812296, -1.6452603299701125, 1.6874626980847942]' OUTTAG=s2_exp_back_le_0 timeout 2400 python3 cmb4.py exp back le 30 > gr2/v4/s2/s2_exp_back_le_0.log 2>&1
