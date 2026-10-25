cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.03847525975408, -1.516406680605399, 0.6308923719165043, 1.4626956399431748, -0.2415476715984184, 4.699789502271772]' OUTTAG=s2_soft_fwd_le_2 timeout 2400 python3 cmb4.py soft fwd le 30 > gr2/v4/s2/s2_soft_fwd_le_2.log 2>&1
