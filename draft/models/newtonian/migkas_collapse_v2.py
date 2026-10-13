"""Chris's exact sheet formula: A/A_all = |s_k kR n + dt (g_gal - g_obs)| / (kR), s_k=+1 expand, -1 collapse.
Point mass at CoE, D=58.2, R=14.3. eps = GM dt/(D^2 kR)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import argrelextrema as ae
D,R=58.2,14.3; ang=np.linspace(0,360,2881); t=np.radians(ang); n=np.stack([-np.cos(t),np.sin(t)],1)
xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1)
dg=(-xp/rp[:,None]**3+xo/D**3)*D**2
def vec(sk,eps): return sk*n+eps*dg
mag=lambda sk,e:np.linalg.norm(vec(sk,e),axis=1); los=lambda sk,e:(vec(sk,e)*n).sum(1)
def ext(y):
    out=[]
    for f,lab in ((np.greater,"max"),(np.less,"min")):
        i=list(ae(np.r_[y[-2],y,y[1]],f)[0]-1); out+= [f"{lab} {y[j]:.3f}@{ang[j]:.0f}°" for j in i if 0<=j<len(y)-1]
    return out
for e in (0.247,1,2):
    print(e,"EXP",ext(mag(1,e)),"| COL",ext(mag(-1,e)),"| COL los range",los(-1,e).min().round(3),los(-1,e).max().round(3))
    print("   shapes differ? max|(exp-1)-(col-1)| =",np.abs(mag(1,e)-mag(-1,e)).max().round(3), " exp-col corr",np.corrcoef(mag(1,e),mag(-1,e))[0,1].round(3))
f,ax=plt.subplots(figsize=(9,5))
ax.axhline(1,color="0.5",lw=1.8,ls="--",label="no mass (= 1)")
for e,lw,a in ((0.247,2.4,1),(1,1.4,.6)):
    ax.plot(ang,mag(1,e),color="#c0392b",lw=lw,alpha=a,label=f"expanding, strength {e}")
    ax.plot(ang,mag(-1,e),color="#2e5fa8",lw=lw,alpha=a,label=f"collapsing, strength {e}")
ax.plot(ang,los(-1,.247),color="#2e5fa8",lw=1.2,ls=":",label="collapsing, signed LOS part (0.247)")
ax.axhline(0,color="k",lw=.6)
for x in (90,270): ax.axvline(x,color="0.7",lw=.8,ls=":")
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle from CoE direction (deg)",fontsize=11)
ax.set_ylabel("A / A_all = |±kR n + dt Δg| / kR",fontsize=11)
ax.set_title("Chris's sheet formula: expanding vs collapsing (point mass, D = 58.2, R = 14.3 Gpc)",fontsize=12)
ax.legend(frameon=False,fontsize=9,ncol=3,loc="upper center",bbox_to_anchor=(0.5,-0.13)); f.savefig("migkas_collapse_v2.png",dpi=130,bbox_inches="tight")
