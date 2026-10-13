import json,glob,numpy as np,sys
for f in sorted(glob.glob(f"raw_{sys.argv[1]}/*.json"),key=lambda f:(len(f),f)):
    d=json.load(open(f)); r=[q for q in d["runs"] if q["viol"]<=1e-3]; L=np.array(sorted(q["left"] for q in r))
    b=min(r,key=lambda q:q["left"]); M=np.array(b["M"]); k=np.array(list(d["kinds"])); X=np.array(b["X"]); dd=np.linalg.norm(X,axis=1)
    print(f"{d['kinds']:14s} best={L[0]:7.2f} med={np.median(L):6.1f} feas={len(L)} Csum={M[k=='C'].sum():.2e} Vsum={M[k=='V'].sum():.2e} d={np.round(dd,1).tolist()} surr={b['left_surroundcomp']:.1f} none={b['left_nocomp']:.1f}")
