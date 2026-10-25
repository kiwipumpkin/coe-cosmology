cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.903522062622244, -2.4044613329882654, 1.6248688618553322, 1.331748530580706, 4.424438513823619, 5.782789733964496]' OUTTAG=s2_soft_fwd_gtc_3 timeout 2400 python3 cmb4.py soft fwd gtc 30 > gr2/v4/s2/s2_soft_fwd_gtc_3.log 2>&1
