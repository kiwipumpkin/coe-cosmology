cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[53.257292435080046, -2.2545477202351125, 1.9257112420064848, 1.5072915221890872, 5.526208890610208, 3.875117980647378]' OUTTAG=s2_soft_fwd_gtc_1 timeout 2400 python3 cmb4.py soft fwd gtc 30 > gr2/v4/s2/s2_soft_fwd_gtc_1.log 2>&1
