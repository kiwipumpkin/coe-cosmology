cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.34203376768829, -1.6063174280887718, 0.15081657671305693, 0.9593118311738891, -0.8806240013661224, 4.2949349457261565]' OUTTAG=s2_soft_fwd_le_4 timeout 2400 python3 cmb4.py soft fwd le 30 > gr2/v4/s2/s2_soft_fwd_le_4.log 2>&1
