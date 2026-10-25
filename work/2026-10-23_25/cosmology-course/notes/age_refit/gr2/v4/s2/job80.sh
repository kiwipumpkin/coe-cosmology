cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.46196038589485333, 1.1296802016354957, 38.02953732712808, 0.09838644263029406, 2.491759082003834, 2.9670946064992827]' OUTTAG=s2_soft_back_le_2 timeout 2400 python3 cmb4.py soft back le 30 > gr2/v4/s2/s2_soft_back_le_2.log 2>&1
