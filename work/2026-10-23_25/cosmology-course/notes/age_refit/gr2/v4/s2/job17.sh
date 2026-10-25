cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[54.94025855460859, -1.7023871960569608, 2.047713766503687, -0.3496077001842979, 3.3099851991770004, 2.1231808939602703]' OUTTAG=s2_eta_fwd_gtc_5 timeout 2400 python3 cmb4.py eta fwd gtc 30 > gr2/v4/s2/s2_eta_fwd_gtc_5.log 2>&1
