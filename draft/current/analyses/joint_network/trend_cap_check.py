"""What leftover is left if the density-trend cap is actually honoured? Refit the strict compromise layouts with the
rule weight raised 10x and 100x (everything else unchanged)."""
import json,numpy as np
from jn_fast import JF,fit
from jn_core import *
from jn_report import score
out={}
for N in (8,24):
    r=[x for x in json.load(open(f"raw/strict_m_N{N}.json")) if x["ok"]]; b=min(r,key=score); p=np.array(b["p"])
    for mult in (10,100):
        m=JF(b["kinds"],"strict",True); p2=fit(m,p,lam=30*mult,stages=(1,1),nfev=600)
        dec=m.decode(p2); ex=JN(b["kinds"],"strict",True).trend_excess(dec); h,a,s=m.split(dec)
        out[f"N{N}_x{mult}"]=dict(left=m.leftover(p2),trend_max_dex=float(ex.max()),soft_max=float(s.max()))
        print(N,mult,out[f"N{N}_x{mult}"],flush=True)
        json.dump(out,open("trend_cap_check.json","w"),indent=1)
