cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/best XOVR='[0.2840633296538337, 1.1914734450372575, 37.934816416812296, -1.6452603299701125, 1.6874626980847942]' OUTTAG=best_exp_back_le_0 python3 cmb4.py exp back le 6 > gr2/v4/best/best_exp_back_le_0.log 2>&1
