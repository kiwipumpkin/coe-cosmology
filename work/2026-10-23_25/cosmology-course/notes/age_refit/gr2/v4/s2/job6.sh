cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.28841612710514175, 1.1914173355277309, 46.18948401305302, 2.9983025505212306, -0.1890260826005753, 3.9515048456497004]' OUTTAG=s2_eta_back_le_0 timeout 2400 python3 cmb4.py eta back le 30 > gr2/v4/s2/s2_eta_back_le_0.log 2>&1
