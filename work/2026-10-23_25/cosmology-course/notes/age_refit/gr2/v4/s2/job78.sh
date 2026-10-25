cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.2840633088410681, 1.1914743801958718, 37.905094529973894, 0.37316625878840987, 2.492527196816506, 3.466288621316759]' OUTTAG=s2_soft_back_le_0 timeout 2400 python3 cmb4.py soft back le 30 > gr2/v4/s2/s2_soft_back_le_0.log 2>&1
