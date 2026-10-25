cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.4663132041589271, 1.1296231569673547, 46.31392681020721, 2.7235227343631148, -0.18979419741324755, 3.452310830832224]' OUTTAG=s2_eta_back_le_1 timeout 2400 python3 cmb4.py eta back le 30 > gr2/v4/s2/s2_eta_back_le_1.log 2>&1
