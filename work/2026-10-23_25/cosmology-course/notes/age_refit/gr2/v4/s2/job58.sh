cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[1.2059426504133859, 0.797925434355151, 48.68429140100523, 3.6030518320690814, -0.32155038074533593, 0.2226642721181067]' OUTTAG=s2_pl_back_le_4 timeout 2400 python3 cmb4.py pl back le 30 > gr2/v4/s2/s2_pl_back_le_4.log 2>&1
