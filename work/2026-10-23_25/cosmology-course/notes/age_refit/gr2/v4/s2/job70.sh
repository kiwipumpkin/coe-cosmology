cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[51.38821517845692, -1.6588320046760074, 2.4411614764867307, 1.8911085036390318, 2.5539123300157973, 0.4047047244351307]' OUTTAG=s2_pl_fwd_le_4 timeout 2400 python3 cmb4.py pl fwd le 30 > gr2/v4/s2/s2_pl_fwd_le_4.log 2>&1
