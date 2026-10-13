"""Which rules stay broken: per-category penalty for a layout (numpy reference code, same order as jn_core.net_pen)."""
import json,numpy as np
from jn_core import *
def cats(kinds,level,mirror,p):
    m=JN(kinds,level,mirror); dec=m.decode(np.array(p)); c=m.isC.sum(); v=m.isV.sum(); n=m.n; out={}
    L=m.L
    if L["cont"]:
        cp=m.cont_pen(dec,L["cont"]); names=["count","void share","mass sizes","void sizes","net contrast"]
        for z,zn in enumerate(["10-14.3 vs 14.3-19","0-14.3 vs 14.3-18"]):
            for i,nm in enumerate(names): out[f"continuity {zn}: {nm}"]=float(cp[5*z+i])
    if L["net"]:
        q=m.net_pen(dec,L["net"]); pv=v*(v-1)//2; seg=[("mass on a void wall",c),("mass between 2 voids",c),("void overlap",pv),("void contact without a mass",pv),
             ("void neighbours 5-7",v),("mass links >=3",c),("mass spacing 1-2 void diameters",c),("void deficit balanced by masses",v),("neighbourhood nets to 0",n)]
        k=0
        for nm,ln in seg: out[nm]=float(q[k:k+ln].max(initial=0)); out[nm+" (n broken)"]=int((q[k:k+ln]>0.1).sum()); k+=ln
        assert k==len(q),(k,len(q))
    if L["trend"]:
        ex=m.trend_excess(dec); out["density trend: max dex above"]=float(ex.max()); out["density trend: n above tol"]=int((ex>L["trend"]["tol_dex"]+0.1).sum())
    return out
if __name__=="__main__":
    import sys,glob
    from jn_report import score
    res={}
    for f in sorted(glob.glob("raw/*.json")):
        r=[x for x in json.load(open(f)) if x["ok"]]
        if not r or r[0]["level"]=="free": continue
        b=min(r,key=score); res[f[4:-5]]=dict(left=b["left"],**cats(b["kinds"],b["level"],b["mirror"],b["p"]))
    json.dump(res,open("rule_breakdown.json","w"),indent=1)
    for k in ["web_m_N24","strict_m_N24","loose_m_N24","cont_m_N24"]:
        print(k,{a:round(b,2) for a,b in res[k].items() if (isinstance(b,float) and b>0.05) or (isinstance(b,int) and b>0)})
