cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[54.90005160445071, -2.658469433941363, 1.048453498329575, 1.8229990502071258, 3.0244215169710253, 2.592527098732862]' OUTTAG=s2_eta_fwd_gtc_3 timeout 2400 python3 cmb4.py eta fwd gtc 30 > gr2/v4/s2/s2_eta_fwd_gtc_3.log 2>&1
