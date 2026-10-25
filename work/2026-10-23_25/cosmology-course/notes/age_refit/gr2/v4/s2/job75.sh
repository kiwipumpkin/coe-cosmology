cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.6507397158661543, 0.8770543985156485, 47.80626148129522, -1.6492590861736671, -1.7351478788925432, 6.273228243158249]' OUTTAG=s2_soft_back_gtc_3 timeout 2400 python3 cmb4.py soft back gtc 30 > gr2/v4/s2/s2_soft_back_gtc_3.log 2>&1
