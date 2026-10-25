cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.3412720444271027, 1.2239497095434138, 46.77814553915232, -1.2339651542126076, -1.3190972792632554, 7.2102806197263565]' OUTTAG=s2_eta_back_gtc_2 timeout 2400 python3 cmb4.py eta back gtc 30 > gr2/v4/s2/s2_eta_back_gtc_2.log 2>&1
