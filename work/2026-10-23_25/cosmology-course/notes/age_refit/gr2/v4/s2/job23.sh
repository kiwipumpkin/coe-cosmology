cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.428364093740484, -2.567067931460146, 0.369459578033418, -3.2711495631382235, -0.7099548373258744, 3.9034172546828954]' OUTTAG=s2_eta_fwd_le_5 timeout 2400 python3 cmb4.py eta fwd le 30 > gr2/v4/s2/s2_eta_fwd_le_5.log 2>&1
