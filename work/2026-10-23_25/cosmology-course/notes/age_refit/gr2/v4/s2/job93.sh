cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.78694471864702, -1.6087623929223294, 0.5828575344470457, -0.02150523812562999, -1.9043283020006854, 4.235134792127967]' OUTTAG=s2_soft_fwd_le_3 timeout 2400 python3 cmb4.py soft fwd le 30 > gr2/v4/s2/s2_soft_fwd_le_3.log 2>&1
