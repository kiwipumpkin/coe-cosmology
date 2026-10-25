cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[54.04230817365533, -2.1915121745538353, 1.2803372461415117, 1.1242104084341902, 0.0012563635922702154, 2.476415786496693]' OUTTAG=s2_pl_fwd_gtc_2 timeout 2400 python3 cmb4.py pl fwd gtc 30 > gr2/v4/s2/s2_pl_fwd_gtc_2.log 2>&1
