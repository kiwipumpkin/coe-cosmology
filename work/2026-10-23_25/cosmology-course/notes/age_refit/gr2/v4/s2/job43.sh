cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.80767432097841, -1.774551747636007, 0.4587154438073757, -0.7986108595344121, -1.336856241628402]' OUTTAG=s2_exp_fwd_le_1 timeout 2400 python3 cmb4.py exp fwd le 30 > gr2/v4/s2/s2_exp_fwd_le_1.log 2>&1
