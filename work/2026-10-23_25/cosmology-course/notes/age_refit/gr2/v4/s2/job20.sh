cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.90863455026319, -1.2078449708017747, 0.24185493980708772, -2.375675886404551, -0.4809518298021589, 4.428866472911355]' OUTTAG=s2_eta_fwd_le_2 timeout 2400 python3 cmb4.py eta fwd le 30 > gr2/v4/s2/s2_eta_fwd_le_2.log 2>&1
