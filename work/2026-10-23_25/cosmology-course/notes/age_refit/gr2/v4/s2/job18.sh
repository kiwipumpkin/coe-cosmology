cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.578978653390905, -1.8047791785168688, 1.1682962171843951, -2.385604412761367, -0.11469554620600086, 3.696139053326365]' OUTTAG=s2_eta_fwd_le_0 timeout 2400 python3 cmb4.py eta fwd le 30 > gr2/v4/s2/s2_eta_fwd_le_0.log 2>&1
