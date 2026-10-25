cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/s2 XOVR='[0.9391558276385914, 0.885745255271686, 47.687027847202856, 3.8865590873179197, 0.8709361071557584, 3.9756674193895867]' OUTTAG=s2_eta_back_le_5 timeout 2400 python3 cmb4.py eta back le 30 > gr2/v4/s2/s2_eta_back_le_5.log 2>&1
