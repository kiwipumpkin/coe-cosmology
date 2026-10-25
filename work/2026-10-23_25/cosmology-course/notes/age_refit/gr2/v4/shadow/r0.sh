cd /workspace/cosmology-course/notes/age_refit
ulimit -v 5000000
XOVR='[54.397677164019115, -2.390068343885659, 1.0872398779500079, 0.25982605975090234, 3.908527530894909, 1.753253085317883]' KFIX=8.460934014852015 OUTTAG=sh_eta_fwd_gtc_1 python3 shadow.py eta fwd gtc > gr2/v4/shadow/sh_eta_fwd_gtc_1.log 2>&1
