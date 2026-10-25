cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.41830022455590155, 1.1442772767511122, 47.62512534676387, 2.837387770537657, -0.508222905294403, 0.3925628118996771]' OUTTAG=s2_pl_back_le_1 timeout 2400 python3 cmb4.py pl back le 30 > gr2/v4/s2/s2_pl_back_le_1.log 2>&1
