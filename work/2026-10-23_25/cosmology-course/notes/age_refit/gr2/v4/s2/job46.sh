cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.713542461027686, -1.507234094002353, 0.7695086042682651, -1.2000777218329604, -0.9115866686327989]' OUTTAG=s2_exp_fwd_le_4 timeout 2400 python3 cmb4.py exp fwd le 30 > gr2/v4/s2/s2_exp_fwd_le_4.log 2>&1
