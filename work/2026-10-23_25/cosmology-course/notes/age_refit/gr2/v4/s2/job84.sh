cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.95534300668561, -1.6132492892633112, 1.343017289870743, 1.1095387203013962, 4.691098327756693, 4.8233607381167936]' OUTTAG=s2_soft_fwd_gtc_0 timeout 2400 python3 cmb4.py soft fwd gtc 30 > gr2/v4/s2/s2_soft_fwd_gtc_0.log 2>&1
