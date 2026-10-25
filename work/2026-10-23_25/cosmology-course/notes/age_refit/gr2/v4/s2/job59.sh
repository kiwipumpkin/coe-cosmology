cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.4090323660960421, 1.3302296539768577, 49.88369967788458, 3.3405361846089363, 0.9237654656171799, 2.1996112829074264]' OUTTAG=s2_pl_back_le_5 timeout 2400 python3 cmb4.py pl back le 30 > gr2/v4/s2/s2_pl_back_le_5.log 2>&1
