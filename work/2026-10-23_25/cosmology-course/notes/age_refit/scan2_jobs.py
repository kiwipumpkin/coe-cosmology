import json,glob,numpy as np,sys
NT=int(sys.argv[1]) if len(sys.argv)>1 else 6
jobs=[]
for f in sorted(glob.glob('gr2/v4/scan1_*.json')):
    d=json.load(open(f)); R=sorted(d['res'],key=lambda r:r['chi2']); ch=[]
    for r in R:
        if len(ch)>=NT: break
        if all(np.max(np.abs(np.array(r['x'])-np.array(c['x'])))>0.1 for c in ch): ch.append(r)
    for i,r in enumerate(ch):
        ot=f"s2_{d['law']}_{d['mode']}_{d['vm']}_{i}"
        jobs.append(f"OUTD=gr2/v4/s2 XOVR='{json.dumps(r['x'])}' OUTTAG={ot} timeout 2400 python3 cmb4.py {d['law']} {d['mode']} {d['vm']} 30 > gr2/v4/s2/{ot}.log 2>&1")
open('gr2/v4/jobs_s2.txt','w').write('\n'.join(jobs)+'\n'); print(len(jobs))
