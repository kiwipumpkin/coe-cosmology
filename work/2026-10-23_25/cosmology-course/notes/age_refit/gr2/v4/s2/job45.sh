cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.9592477433433, -2.1941224842230547, 1.2239833877707842, -0.3486566515478189, 0.0012255922780697315]' OUTTAG=s2_exp_fwd_le_3 timeout 2400 python3 cmb4.py exp fwd le 30 > gr2/v4/s2/s2_exp_fwd_le_3.log 2>&1
