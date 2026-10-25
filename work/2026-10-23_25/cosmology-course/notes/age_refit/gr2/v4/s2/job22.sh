cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[48.796425722539425, -2.25005279481411, 1.4174272594973445, 2.2999815866351128, -2.934197813272476, 2.003122834302485]' OUTTAG=s2_eta_fwd_le_4 timeout 2400 python3 cmb4.py eta fwd le 30 > gr2/v4/s2/s2_eta_fwd_le_4.log 2>&1
