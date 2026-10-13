"""Collapsing A/A_all = |-kR n + dt(g_gal-g_obs)|/kR for point masses = f x sheet mass (4.44e68 kg, strength 0.247)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import argrelextrema as ae
D,R=58.2,14.3; ang=np.linspace(0,360,2881); t=np.radians(ang); n=np.stack([-np.cos(t),np.sin(t)],1)
xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1); dg=(-xp/rp[:,None]**3+xo/D**3)*D**2
F=[0,0.5,1,2,4,8]; cm=plt.cm.viridis(np.linspace(0.05,0.9,len(F)))
def ext(y):
    o=[]
    for f,l in ((np.greater,"peak"),(np.less,"dip")):
        o+=[f"{l} {y[j]:+.2f}@{ang[j]:.0f}°" for j in ae(np.r_[y[-2],y,y[1]],f)[0]-1 if 0<=j<len(y)-1 and ang[j]<=180]
    return o
for kind in ("mag","signed"):
    f,ax=plt.subplots(figsize=(9,5))
    ax.axhline(1 if kind=="mag" else -1,color="0.5",lw=1.6,ls="--",label="no mass")
    if kind=="signed": ax.axhline(0,color="k",lw=1)
    for fm,c in zip(F[1:],cm[1:]):
        v=-n+0.247*fm*dg; y=np.linalg.norm(v,axis=1) if kind=="mag" else (v*n).sum(1)
        ax.plot(ang,y,color=c,lw=2,label=f"{fm}× : M = {fm*4.44:.3g}×10⁶⁸ kg, strength {0.247*fm:.3g}")
        print(kind,fm,ext(y),"receding frac",round(float(np.mean(y>0)),2) if kind=="signed" else "")
    for x in (90,270): ax.axvline(x,color="0.75",lw=.8,ls=":")
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle from CoE direction (deg)",fontsize=11)
    ax.set_ylabel("A / A_all = |−kR n + dt Δg| / kR" if kind=="mag" else "signed LOS part / kR  (+ = receding)",fontsize=11)
    ax.set_title("Collapsing (V = −k r), point mass at CoE, D = 58.2, R = 14.3 Gpc"+("" if kind=="mag" else ": signed"),fontsize=12)
    ax.legend(frameon=False,fontsize=9,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.13))
    f.savefig("migkas_collapse_masses.png" if kind=="mag" else "migkas_collapse_masses_signed.png",dpi=130,bbox_inches="tight"); plt.close(f)
