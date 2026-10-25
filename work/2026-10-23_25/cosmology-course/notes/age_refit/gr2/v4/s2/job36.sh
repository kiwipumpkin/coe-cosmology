cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[56.55353603461143, -1.7052106735557893, 0.45759675316127935, 0.13688216939813627, 0.024667291054779862]' OUTTAG=s2_exp_fwd_gtc_0 timeout 2400 python3 cmb4.py exp fwd gtc 30 > gr2/v4/s2/s2_exp_fwd_gtc_0.log 2>&1
