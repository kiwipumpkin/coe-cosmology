cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[55.05059066652393, -2.198474778443241, 1.1217264337531512, 0.29470728292308046, -0.08333786005603938, 2.849134955159718]' OUTTAG=s2_pl_fwd_gtc_0 timeout 2400 python3 cmb4.py pl fwd gtc 30 > gr2/v4/s2/s2_pl_fwd_gtc_0.log 2>&1
