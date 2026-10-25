cd /workspace/cosmology-course/notes/age_refit
ulimit -v 5000000
[ -f gr2/v4/shadow/refit_lM53.83.json ] || python3 refitM.py 53.826248592590545 > gr2/v4/mscan/refit_lM53.826.log 2>&1
X=$(python3 -c "import json;print(json.dumps(json.load(open('gr2/v4/shadow/refit_lM53.83.json'))['x']))")
[ -f gr2/v4/mscan/lM53.826.json ] || OUTD=gr2/v4/mscan XOVR="$X" OUTTAG=lM53.826 python3 cmb4.py eta fwd gtc 6 > gr2/v4/mscan/lM53.826.log 2>&1
K=$(python3 -c "import json;print(json.load(open('gr2/v4/mscan/lM53.826.json'))['K_local_today_K'])")
XOVR="$X" KFIX=$K OUTTAG=shs_lM53.826 python3 shadow.py eta fwd gtc > gr2/v4/mscan/shs_lM53.826.log 2>&1
cp gr2/v4/shadow/shs_lM53.826.json gr2/v4/mscan/ 2>/dev/null
