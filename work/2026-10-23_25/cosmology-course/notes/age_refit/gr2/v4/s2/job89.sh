cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.665516529232264, -2.232803001999855, -0.09803781099617481, 3.8541608676314354, -2.6800486147403717, 3.1395473398268225]' OUTTAG=s2_soft_fwd_gtc_5 timeout 2400 python3 cmb4.py soft fwd gtc 30 > gr2/v4/s2/s2_soft_fwd_gtc_5.log 2>&1
