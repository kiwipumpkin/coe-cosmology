"""Signed version of Chris's sheet quantity: S(th) = [ s_k kR n + dt (g_P - g_obs) ] . n / (kR),
s_k=+1 expanding (V=kr), -1 collapsing (V=-kr). eps = GM dt/(D^2 kR). D=58.2, R=14.3, point mass at CoE.
Positive = receding. Also free-fall collapse v=-sqrt(2GM/r) (no dt step) for comparison."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D,R=58.2,14.3
ang=np.linspace(0,360,2881); t=np.radians(ang); n=np.stack([-np.cos(t),np.sin(t)],1)
xo=np.array([D,0.])
def at(Rs):
    xp=xo+Rs*n; rp=np.linalg.norm(xp,axis=1)
    dg=(-xp/rp[:,None]**3+xo/D**3)*D**2; return (dg*n).sum(1),xp,rp
ta,xp,rp=at(R)
S=lambda sk,eps:sk+eps*ta
print("tilde a: max",ta.max(),"at",ang[ta.argmax()],"min",ta.min(),"at",ang[ta.argmin()])
print("collapse: recession appears where eps*ta>1; first at eps=",1/ta.max())
for eps in [0.247,1,1.5,2,3,5]:
    y=S(-1,eps); pos=np.mean(y>0); i=np.where(np.diff(np.sign(y))!=0)[0]
    ye=S(1,eps); print(f"eps={eps}: collapse receding fraction of circle {pos:.2f}, flips {np.round(ang[i],1)}; expanding min {ye.min():.3f}")
print("expanding band goes negative for eps >",1/-ta.min())
# free fall collapse, inner shells faster: v=-sqrt(2GM/r)
v=lambda r:-np.sqrt(1/r); vr=(v(rp)[:,None]*xp/rp[:,None]-v(D)*xo/D); vl=(vr*n).sum(1)
i=np.where(np.diff(np.sign(vl))!=0)[0]; print("free-fall collapse: receding fraction",np.mean(vl>0),"flips",np.round(ang[i],1))
def fig(curves,title,fn):
    f,ax=plt.subplots(figsize=(8.5,4.6))
    for y,c,ls,l in curves: ax.plot(ang,y,color=c,lw=2.2,ls=ls,label=l)
    ax.axhline(0,color="k",lw=0.8)
    lo,hi=ax.get_ylim(); ax.axhspan(0,hi,color="#fbe3e3",zorder=0); ax.axhspan(lo,0,color="#e3ebfb",zorder=0); ax.set_ylim(lo,hi)
    for a in (90,270): ax.axvline(a,color="0.6",lw=0.8,ls=":")
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle from CoE direction (deg)",fontsize=11)
    ax.set_ylabel("signed v_los with gravity / kR   (+ = receding)",fontsize=11); ax.set_title(title,fontsize=12)
    ax.legend(frameon=False,fontsize=10); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
fig([(S(1,.247),"#c0392b","-","expanding, V = +k r"),(S(-1,.247),"#2e5fa8","-","collapsing, V = −k r"),(0*ang+1,"0.5","--","no mass (±1)"),(0*ang-1,"0.5","--",None)],
    "Chris's sheet geometry, point mass, GM·dt/D²kR = 0.247","migkas_collapse.png")
fig([(S(-1,e),c,"-",f"collapsing, GM·dt/D²kR = {e}") for e,c in ((1.32,"#7f8fa6"),(2,"#9b59b6"),(3,"#c0392b"))],
    "Collapse with strong gravity step: recession except near 90°/270°","migkas_collapse_condition.png")
