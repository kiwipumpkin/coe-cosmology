cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[55.21414293856677, -2.163099780611751, -0.1511392217167935, 1.046424068736501, -1.3357905446642095]' OUTTAG=s2_exp_fwd_gtc_3 timeout 2400 python3 cmb4.py exp fwd gtc 30 > gr2/v4/s2/s2_exp_fwd_gtc_3.log 2>&1
