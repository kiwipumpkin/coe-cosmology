cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/best XOVR='[54.397677164019115, -2.390068343885659, 1.0872398779500079, 0.25982605975090234, 3.908527530894909, 1.753253085317883]' OUTTAG=best_eta_fwd_gtc_1 python3 cmb4.py eta fwd gtc 6 > gr2/v4/best/best_eta_fwd_gtc_1.log 2>&1
