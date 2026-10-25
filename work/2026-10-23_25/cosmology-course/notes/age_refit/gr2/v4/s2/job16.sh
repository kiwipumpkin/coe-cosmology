cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.89109579050637, -1.9237097232122409, 2.2731736960593416, 1.037747218571293, 5.106618279756341, 1.7339161576357216]' OUTTAG=s2_eta_fwd_gtc_4 timeout 2400 python3 cmb4.py eta fwd gtc 30 > gr2/v4/s2/s2_eta_fwd_gtc_4.log 2>&1
