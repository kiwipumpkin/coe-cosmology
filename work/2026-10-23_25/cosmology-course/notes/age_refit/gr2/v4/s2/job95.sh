cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[51.644607179413754, -1.5648894413924093, 0.2810751210612372, 1.7787585643953008, -0.8161502457084495, 3.976365257599712]' OUTTAG=s2_soft_fwd_le_5 timeout 2400 python3 cmb4.py soft fwd le 30 > gr2/v4/s2/s2_soft_fwd_le_5.log 2>&1
