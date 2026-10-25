cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[54.6481372709859, -2.3181620804813146, 0.3858123309976642, -0.14761040498324451, -0.17959413571732918, 2.9600579976805985]' OUTTAG=s2_pl_fwd_gtc_4 timeout 2400 python3 cmb4.py pl fwd gtc 30 > gr2/v4/s2/s2_pl_fwd_gtc_4.log 2>&1
