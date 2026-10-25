cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[52.72684633112421, -0.766320741789345, 0.9347313476836985, -3.6908365958378084, -0.3564407648406063, 4.466364893015739]' OUTTAG=s2_eta_fwd_le_3 timeout 2400 python3 cmb4.py eta fwd le 30 > gr2/v4/s2/s2_eta_fwd_le_3.log 2>&1
