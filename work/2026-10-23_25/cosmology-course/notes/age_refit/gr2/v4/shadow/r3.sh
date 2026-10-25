cd /workspace/cosmology-course/notes/age_refit
ulimit -v 5000000
python3 refitM.py 52.397677164019115 > gr2/v4/shadow/refit_-2.log 2>&1
X=$(python3 -c "import json;print(json.dumps(json.load(open('gr2/v4/shadow/refit_lM52.40.json'))['x']))")
XOVR="$X" KFIX=8.460934014852015 OUTTAG=shM_-2 python3 shadow.py eta fwd gtc > gr2/v4/shadow/shM_-2.log 2>&1
