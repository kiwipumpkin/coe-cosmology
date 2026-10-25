cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.6296881561995398, 1.2326405662994513, 46.65891190505995, 4.30185301927898, 1.2652705398247424, 4.912719795957904]' OUTTAG=s2_eta_back_le_3 timeout 2400 python3 cmb4.py eta back le 30 > gr2/v4/s2/s2_eta_back_le_3.log 2>&1
