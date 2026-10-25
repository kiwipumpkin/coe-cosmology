cd /workspace/cosmology-course/notes/age_refit
ulimit -v 5000000
XOVR='[52.95534300668561, -1.6132492892633112, 1.343017289870743, 1.1095387203013962, 4.691098327756693, 4.8233607381167936]' KFIX=19.87736245404969 OUTTAG=sh_soft_fwd_gtc_0 python3 shadow.py soft fwd gtc > gr2/v4/shadow/sh_soft_fwd_gtc_0.log 2>&1
