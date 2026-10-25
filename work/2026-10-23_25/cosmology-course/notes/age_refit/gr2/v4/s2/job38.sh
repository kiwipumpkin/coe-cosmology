cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[59.229344855993986, -1.5402775928378105, 3.658739998936653, -3.667526535689831, 2.4425913393497467]' OUTTAG=s2_exp_fwd_gtc_2 timeout 2400 python3 cmb4.py exp fwd gtc 30 > gr2/v4/s2/s2_exp_fwd_gtc_2.log 2>&1
