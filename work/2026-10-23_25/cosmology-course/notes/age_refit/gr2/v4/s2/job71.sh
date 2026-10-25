cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[49.46392992867233, -1.5282573256738416, 2.0014913412196518, 4.050498470737864, 0.24595774483611077, 1.6295178001712958]' OUTTAG=s2_pl_fwd_le_5 timeout 2400 python3 cmb4.py pl fwd le 30 > gr2/v4/s2/s2_pl_fwd_le_5.log 2>&1
