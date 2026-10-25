cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.24040314750211622, 1.2060714553114884, 47.500682549609685, 3.1121675866957728, -0.507454790481731, 0.8917568267171534]' OUTTAG=s2_pl_back_le_0 timeout 2400 python3 cmb4.py pl back le 30 > gr2/v4/s2/s2_pl_back_le_0.log 2>&1
