cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[54.80013055955715, -2.2703810418475854, 1.8231539807054948, 0.7021437476572272, 4.0047838065561985, 1.6423300427970022]' OUTTAG=s2_eta_fwd_gtc_0 timeout 2400 python3 cmb4.py eta fwd gtc 30 > gr2/v4/s2/s2_eta_fwd_gtc_0.log 2>&1
