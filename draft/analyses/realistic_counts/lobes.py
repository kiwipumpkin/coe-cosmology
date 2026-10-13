import json, numpy as np
from counts import CL, VO, SCL, SVO, RHO, GPC
z=np.load("field_inside_mirror.npz"); d=z["arr_0"]; P=z["arr_1"]; H=0.3
m=np.full(len(d),RHO*(H*GPC)**3)          # ~mean density (exp-lock varies <2%)
lon=np.degrees(np.arctan2(P[:,1],P[:,0]))%360; r=np.linalg.norm(P,axis=1)
cent={"+":[15,128,195,263],"-":[75,162,228,315]}
out=[]
for sg,cs in cent.items():
    sel=(d>0.005) if sg=="+" else (d<-0.005)
    for c in cs:
        dl=np.abs((lon-c+180)%360-180)
        other=np.min([np.abs((lon-o+180)%360-180) for o in cs],0)
        k=sel&(dl==other)&(dl<45)
        M=float(np.abs(m[k]*d[k]).sum()); V=float(k.sum()*H**3)
        w=np.abs(d[k]); Pc=(P[k]*w[:,None]).sum(0)/w.sum()
        out.append(dict(sign=sg,lon=c,M=M,V_gpc3=V,mean_delta=float(d[k].mean()),peak=float(np.abs(d[k]).max()),
            centroid=Pc.tolist(),dist=float(np.linalg.norm(Pc)),
            n_clusters=M/CL["excess"] if sg=="+" else 0,n_superclusters=M/SCL["excess"] if sg=="+" else 0,
            n_voids=M/VO["deficit"] if sg=="-" else 0,n_supervoids=M/SVO["deficit"] if sg=="-" else 0,
            void_fill_extra=abs(float(d[k].mean()))/0.4 if sg=="-" else 0))
for o in out: print({k:(round(v,3) if isinstance(v,float) and abs(v)<100 else (f"{v:.2e}" if isinstance(v,float) else v)) for k,v in o.items() if k!="centroid"})
json.dump(out,open("lobes_inside.json","w"),indent=1)
