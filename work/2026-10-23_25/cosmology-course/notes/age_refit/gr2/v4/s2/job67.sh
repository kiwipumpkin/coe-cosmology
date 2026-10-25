cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[49.602668542248026, -1.858724925117933, 1.6770377042157532, 3.8902730244018926, 1.6959442836858818, 1.112200581465284]' OUTTAG=s2_pl_fwd_le_1 timeout 2400 python3 cmb4.py pl fwd le 30 > gr2/v4/s2/s2_pl_fwd_le_1.log 2>&1
