cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.5990406328011497, 1.20908205781306, 38.403066791122015, -0.44406288457750387, 2.084128333603524]' OUTTAG=s2_exp_back_le_2 timeout 2400 python3 cmb4.py exp back le 30 > gr2/v4/s2/s2_exp_back_le_2.log 2>&1
