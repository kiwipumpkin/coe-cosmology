cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[1.5332704726274114e-08, 1.1827264787716933, 46.308717647145386, -2.537515622970356, -2.7733939016885727, 6.249065669418153]' OUTTAG=s2_eta_back_gtc_0 timeout 2400 python3 cmb4.py eta back gtc 30 > gr2/v4/s2/s2_eta_back_gtc_0.log 2>&1
