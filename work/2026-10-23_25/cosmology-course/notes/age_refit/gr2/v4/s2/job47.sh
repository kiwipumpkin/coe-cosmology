cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.413575760652606, -1.4280605086023677, 1.3396683617780927, -1.730633378012746, -1.0141722423134978]' OUTTAG=s2_exp_fwd_le_5 timeout 2400 python3 cmb4.py exp fwd le 30 > gr2/v4/s2/s2_exp_fwd_le_5.log 2>&1
