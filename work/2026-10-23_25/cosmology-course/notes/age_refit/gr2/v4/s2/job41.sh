cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[57.039542674509256, -1.8679999904037095, 0.301331485058991, -1.2585416499128268, 0.8999261698595231]' OUTTAG=s2_exp_fwd_gtc_5 timeout 2400 python3 cmb4.py exp fwd gtc 30 > gr2/v4/s2/s2_exp_fwd_gtc_5.log 2>&1
