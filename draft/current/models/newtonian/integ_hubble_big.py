import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("integ_hik.py").read().split("def calc")[0])     # field(sg,m,kk) RK4, k, T, D, ang, n, w
K=3*k; Ms=[0,1,4,10]; F={m:field(1,m,K)[0] for m in Ms}
xo=np.array([D,0.]); KMSMPC=299792.458/1000
def vlos(m,Rv,nn):
    xp=xo+Rv*nn; rp=np.linalg.norm(xp,axis=-1); vf=F[m]
    vr=vf(rp)[...,None]*xp/rp[...,None]-vf(np.array(D))*xo/D; return (vr*nn).sum(-1)/c
def Rmax(m):
    for Rv in np.arange(14.3,40,0.5):
        if F[m](np.linspace(D-Rv,D+Rv,300)).min()<=0: return Rv-0.5
    return 40
Rlim={m:min(Rmax(m),25) for m in Ms}; print("R limit (outward filter, capped 25):",Rlim)

plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":15,"xtick.labelsize":16,"ytick.labelsize":16})
Rs=np.linspace(0.5,25,200); dirs=[0,45,90,135,180]; cols=["#b2182b","#ef8a62","#4d4d4d","#67a9cf","#2166ac"]; mean=lambda y:(w*y).sum()/w.sum()
for m in Ms:
    f,ax=plt.subplots(figsize=(8,6)); R_=Rs[Rs<=Rlim[m]]
    for a_,cl in zip(dirs,cols):
        nn=np.array([np.cos(np.radians(a_)),np.sin(np.radians(a_))]); ax.plot(R_,[vlos(m,r,nn) for r in R_],color=cl,lw=3.5,label="%d°"%a_)
    ax.axhline(0,color="k",lw=1); ax.set_xlabel("galaxy distance R (Gpc)"); ax.set_ylabel("v_LOS / c"); ax.set_ylim(-0.01,0.47)
    ax.set_title("Expanding, "+("no mass" if m==0 else "M = %g×"%m)); ax.legend(frameon=False,title="direction",loc="upper left")
    f.tight_layout(); f.savefig("integ_hubble_big_diagram_M%g.png"%m,dpi=110); plt.close(f)
h=0.05
def slope(m,Rv): return (vlos(m,Rv+h,n)-vlos(m,Rv-h,n))/(2*h)
cm=plt.cm.viridis(np.linspace(.1,.9,3))
for Rv,tag in ((1.0,"R1"),(14.3,"R14")):
    f,ax=plt.subplots(figsize=(8,6)); j=0
    for m in Ms:
        y=slope(m,Rv)*KMSMPC
        if m==0: ax.plot(ang,y,"k--",lw=2.5,label="no mass")
        else: ax.plot(ang,y,color=cm[j],lw=3.5,label="%g×"%m); j+=1
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,90)); ax.set_ylim(1,7.3); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_ylabel("dv/dR (km/s/Mpc)")
    ax.set_title("Hubble rate vs angle, R = %g Gpc"%Rv); ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18),columnspacing=1)
    f.tight_layout(); f.savefig("integ_hubble_big_slope_%s.png"%tag,dpi=110); plt.close(f)
