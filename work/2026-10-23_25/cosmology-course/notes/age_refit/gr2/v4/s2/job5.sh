cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[8.482230920344591, 0.8127885432913899, 46.33045865222812, -2.4324854984879494, 2.5802564918994904, 4.713840781897306]' OUTTAG=s2_eta_back_gtc_5 timeout 2400 python3 cmb4.py eta back gtc 30 > gr2/v4/s2/s2_eta_back_gtc_5.log 2>&1
