cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.5894993661021, -2.11656043984841, 0.5773063705588701, -0.775210360048979, -0.22088804404424978]' OUTTAG=s2_exp_fwd_le_0 timeout 2400 python3 cmb4.py exp fwd le 30 > gr2/v4/s2/s2_exp_fwd_le_0.log 2>&1
