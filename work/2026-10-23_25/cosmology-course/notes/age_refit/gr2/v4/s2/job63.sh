cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[54.91871096990144, -1.7409903524663686, 2.271207681863684, -1.1156586590112072, -0.18853740288394372, 2.3026253018298686]' OUTTAG=s2_pl_fwd_gtc_3 timeout 2400 python3 cmb4.py pl fwd gtc 30 > gr2/v4/s2/s2_pl_fwd_gtc_3.log 2>&1
