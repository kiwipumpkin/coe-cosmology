cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[55.72414488769209, -2.2018739111941694, -0.0514823226191814, -0.8092010323778611, -1.2400510058410732, 3.491907795495683]' OUTTAG=s2_pl_fwd_gtc_5 timeout 2400 python3 cmb4.py pl fwd gtc 30 > gr2/v4/s2/s2_pl_fwd_gtc_5.log 2>&1
