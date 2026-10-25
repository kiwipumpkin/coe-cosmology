cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.60976164261719, -1.634868540855054, 1.9239277972019022, -2.439414079767502, 0.8252522880324091, 3.3633527074719636]' OUTTAG=s2_eta_fwd_le_1 timeout 2400 python3 cmb4.py eta fwd le 30 > gr2/v4/s2/s2_eta_fwd_le_1.log 2>&1
