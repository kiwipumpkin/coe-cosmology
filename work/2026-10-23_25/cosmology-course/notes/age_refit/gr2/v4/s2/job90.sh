cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.64681672359558, -2.176756238712954, 0.35829774864885366, 1.030246209715895, -1.2095296946214873, 3.7542839409646995]' OUTTAG=s2_soft_fwd_le_0 timeout 2400 python3 cmb4.py soft fwd le 30 > gr2/v4/s2/s2_soft_fwd_le_0.log 2>&1
