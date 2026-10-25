cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.564744838247095, 0.8845085931214212, 37.62109427293127, -2.431585877814284, 0.4207454969859945]' OUTTAG=s2_exp_back_le_4 timeout 2400 python3 cmb4.py exp back le 30 > gr2/v4/s2/s2_exp_back_le_4.log 2>&1
