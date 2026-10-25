cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.994124495849974, -1.8535005944300877, 3.1966409568447425, 1.4702711359780904, 3.378079079256687, 1.4182962717142433]' OUTTAG=s2_eta_fwd_gtc_2 timeout 2400 python3 cmb4.py eta fwd gtc 30 > gr2/v4/s2/s2_eta_fwd_gtc_2.log 2>&1
