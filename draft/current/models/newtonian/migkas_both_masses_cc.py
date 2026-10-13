"""Expanding/collapsing A/A_all = |±kR n + dt(g_gal-g_obs)|/kR and signed LOS, masses f x sheet (4.44e68 kg, strength 0.247)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import argrelextrema as ae
D,R=58.2,14.3; ang=np.linspace(0,360,2881); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1)
xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1); dg=(-xp/rp[:,None]**3+xo/D**3)*D**2
F=[0,0.5,1,2,4,8]; cm=plt.cm.viridis(np.linspace(0.05,0.9,len(F)))
def Y(sk,fm,kind):
    v=sk*n+0.247*fm*dg; return np.linalg.norm(v,axis=1) if kind=="mag" else (v*n).sum(1)
def ext(y):
    o=[]
    for f,l in ((np.greater,"pk"),(np.less,"dip")):
        o+=[f"{l} {y[j]:+.2f}@{ang[j]:.0f}" for j in ae(np.r_[y[-2],y,y[1]],f)[0]-1 if 0<=j<len(y)-1 and ang[j]<=180]
    i=np.where(np.diff(np.sign(y))!=0)[0]; return o+([f"sign change @{[round(float(ang[k]),1) for k in i]}"] if len(i) else [])
for kind in ("mag","signed"):
    allv=np.concatenate([Y(s,f,kind) for s in (1,-1) for f in F]); lo,hi=allv.min(),allv.max(); pad=.05*(hi-lo)
    for sk,name in ((1,"exp"),(-1,"col")):
        f,ax=plt.subplots(figsize=(9,5))
        if kind=="signed": ax.axhline(0,color="k",lw=1.2)
        for fm,c in zip(F,cm):
            y=Y(sk,fm,kind); print(name,kind,fm,ext(y) if fm else "flat")
            ax.plot(ang,y,color=c,lw=2,ls="--" if fm==0 else "-",label="0× : no mass" if fm==0 else f"{fm}× : M = {fm*4.44:.3g}×10⁶⁸ kg, strength {0.247*fm:.3g}")
        for x in (90,270): ax.axvline(x,color="0.75",lw=.8,ls=":")
        ax.set_ylim(lo-pad,hi+pad); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11)
        ax.set_ylabel("A / A_all = |±kR n + dt Δg| / kR" if kind=="mag" else "signed LOS part / kR  (+ = receding)",fontsize=11)
        ax.set_title(f"{'Expanding (V = +k r)' if sk>0 else 'Collapsing (V = −k r)'}, point mass at CoE"+(", signed" if kind=="signed" else ""),fontsize=12)
        ax.legend(frameon=False,fontsize=9,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.13))
        f.savefig(f"migkas_{name}_masses{'_signed' if kind=='signed' else ''}_cc.png",dpi=130,bbox_inches="tight"); plt.close(f)
