cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/best XOVR='[49.602668542248026, -1.858724925117933, 1.6770377042157532, 3.8902730244018926, 1.6959442836858818, 1.112200581465284]' OUTTAG=best_pl_fwd_le_1 python3 cmb4.py pl fwd le 6 > gr2/v4/best/best_pl_fwd_le_1.log 2>&1
