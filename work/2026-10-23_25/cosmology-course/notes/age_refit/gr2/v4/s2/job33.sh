cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.9362550747711404, 1.1492884188162802, 36.73440875065809, -0.5052917400905401, 2.896239674998008]' OUTTAG=s2_exp_back_le_3 timeout 2400 python3 cmb4.py exp back le 30 > gr2/v4/s2/s2_exp_back_le_3.log 2>&1
