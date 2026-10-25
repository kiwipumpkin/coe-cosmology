cd /workspace/cosmology-course/notes/age_refit
ulimit -v 5000000
XOVR='[49.602668542248026, -1.858724925117933, 1.6770377042157532, 3.8902730244018926, 1.6959442836858818, 1.112200581465284]' KFIX=59.272380099278244 OUTTAG=sh_pl_fwd_le_1 python3 shadow.py pl fwd le > gr2/v4/shadow/sh_pl_fwd_le_1.log 2>&1
