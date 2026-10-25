cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[56.40435494152215, -2.0502941328362203, 0.075682303305594, 0.5182387141907885, -0.6404019630604441]' OUTTAG=s2_exp_fwd_gtc_4 timeout 2400 python3 cmb4.py exp fwd gtc 30 > gr2/v4/s2/s2_exp_fwd_gtc_4.log 2>&1
