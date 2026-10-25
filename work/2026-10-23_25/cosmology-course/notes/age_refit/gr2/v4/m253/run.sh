cd /workspace/cosmology-course/notes/age_refit
ulimit -v 6000000
OUTD=gr2/v4/m253 XOVR='[53.397677164019115, -1.4354674916437993, 1.1049963828041158, 0.2589794506851767, 3.9159418340317687, 1.7552286729686553]' OUTTAG=m253_full python3 cmb4.py eta fwd gtc 6 > gr2/v4/m253/full.log 2>&1
K=$(python3 -c "import json;print(json.load(open('gr2/v4/m253/m253_full.json'))['K_local_today_K'])")
XOVR='[53.397677164019115, -1.4354674916437993, 1.1049963828041158, 0.2589794506851767, 3.9159418340317687, 1.7552286729686553]' KFIX=$K OUTTAG=m253_ring python3 ring.py eta fwd gtc > gr2/v4/m253/ring.log 2>&1
