"""Tube-width convergence: same seeds, cells kept within 3s / 6s (production) / 12s of the CMB ring.
Reports power and projection share of the 6s pattern relative to the 12s pattern, for foam and unstructured."""
import json, sys, numpy as np
from web_core import *
out={}
for s in (0.5,0.25):
    rows=[]
    for i in range(int(sys.argv[1]) if len(sys.argv)>1 else 6):
        ws={}
        for T in (3*s,6*s,12*s):
            X,M,rad,info=foam(s,np.random.default_rng([777,i]),gen=tube_gen(12*s),sel=lambda S:ring_dist(S)<T)
            ws[T]=pattern(X,M)
            if T==12*s: Mb,radb=M,rad
        Xu=unstructured_tube(Mb,radb,np.random.default_rng([778,i]),12*s); d=ring_dist(Xu); wu={}
        for T in (3*s,6*s,12*s):
            k=d<T; wu[T]=pattern(Xu[k],Mb[k])
        ref=ws[12*s]; refu=wu[12*s]
        r=dict(foam_power_3s=float(np.sum(ws[3*s]**2)/np.sum(ref**2)),foam_power_6s=float(np.sum(ws[6*s]**2)/np.sum(ref**2)),
               foam_proj_6s=float(np.dot(ws[6*s],ref)/np.sum(ref**2)),
               unstr_power_3s=float(np.sum(wu[3*s]**2)/np.sum(refu**2)),unstr_power_6s=float(np.sum(wu[6*s]**2)/np.sum(refu**2)),
               unstr_proj_6s=float(np.dot(wu[6*s],refu)/np.sum(refu**2)),foam_rms_12s=rms(ref)/FLOOR,unstr_rms_12s=rms(refu)/FLOOR)
        rows.append(r); print(s,i,{k:round(v,3) for k,v in r.items()},flush=True)
    out[str(s)]=dict(rows=rows,median={k:float(np.median([r[k] for r in rows])) for k in rows[0]})
json.dump(out,open("check_tube.json","w"),indent=1)
