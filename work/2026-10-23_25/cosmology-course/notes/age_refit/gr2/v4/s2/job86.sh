cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.966251795687015, -2.200760874357346, 1.795456273936991, 0.8266809455445474, 6.129748051994314, 4.571617359999113]' OUTTAG=s2_soft_fwd_gtc_2 timeout 2400 python3 cmb4.py soft fwd gtc 30 > gr2/v4/s2/s2_soft_fwd_gtc_2.log 2>&1
