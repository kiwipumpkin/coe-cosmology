cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.28909541257221, -1.6469508645022977, 0.4850518501649621, -1.1581624196681077, -0.024414624377959493]' OUTTAG=s2_exp_fwd_le_2 timeout 2400 python3 cmb4.py exp fwd le 30 > gr2/v4/s2/s2_exp_fwd_le_2.log 2>&1
