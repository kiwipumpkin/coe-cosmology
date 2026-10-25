cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[1.5332704726274114e-08, 1.1827264787716933, 46.308717647145386, -2.537515622970356, -2.795110068648877, 6.249065669418362]' OUTTAG=s2_soft_back_gtc_0 timeout 2400 python3 cmb4.py soft back gtc 30 > gr2/v4/s2/s2_soft_back_gtc_0.log 2>&1
