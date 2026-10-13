"""Point mass 7.5e52 kg at CoE only. Collapse: (a) kappa=1.05, m=-4 stall family (r_max = kappa D (r/D)^m, shells past stall at rest);
(b) far-fall: r_max = inf. X = D (20-80 Gpc), Y = s (1-30 Gpc); cells with s >= 0.9 D masked (sightline passes too near the CoE).
Speeds scale as sqrt(M); fractions are mass-independent."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm,SymLogNorm
G=6.674e-11; GPC=3.0857e25; MP=7.52e52
th=np.radians(np.linspace(0,180,721)); ct,st=np.cos(th),np.sin(th); w=st/st.sum()
def stats(D,s,far):
    def v(r):
        r=np.asarray(r,float); rmax=np.inf if far else 1.05*D*(r/D)**-4
        return -np.sqrt(np.maximum(2*G*MP*(1/r-1/rmax)/GPC,0))
    x=D+s*ct; y=s*st; r=np.hypot(x,y); vl=v(r)*(x*ct+y*st)/r-v(np.array(D))*ct
    return (w*(vl>0)).sum(),(w*vl).sum()/1e3
DD=np.linspace(20,80,61); S=np.linspace(1,30,59)
for far,suf,lab in ((False,"","stall family κ = 1.05, m = −4"),(True,"_farfall","all shells fall from far out (r_max = ∞)")):
    R=np.full((len(S),len(DD)),np.nan); V=R.copy()
    for i,s in enumerate(S):
        for j,D in enumerate(DD):
            if s<0.9*D: R[i,j],V[i,j]=stats(D,s,far)
    print(lab,"| frac %.3f-%.3f | mean v %.0f to %.0f km/s"%(np.nanmin(R),np.nanmax(R),np.nanmin(V),np.nanmax(V)))
    for D in (20,40,58.2,80):
        j=np.argmin(abs(DD-D)); r=R[:,j]; v=V[:,j]; ok=~np.isnan(r)
        c5=S[ok][np.where(np.diff(np.sign(r[ok]-.5)))[0]]; c0=S[ok][np.where(np.diff(np.sign(v[ok])))[0]]
        print(f"   D={DD[j]:.0f}: frac s=1 {r[0]:.2f}, s=max {r[ok][-1]:.2f}; mean {v[0]:.0f} -> {v[ok][-1]:.0f}; 50% at s={np.round(c5,1)}; zero at s={np.round(c0,1)}")
    for Z,norm,cl,lev,fl,tt in ((R*100,TwoSlopeNorm(50,0,100),"sky fraction receding (%)",50,f"heat_pm_frac_D_s{suf}.png",f"Collapse, point mass only — {lab}\nreceding sky fraction (independent of mass)"),
                                (V,SymLogNorm(linthresh=500,vmin=-np.nanmax(abs(V)),vmax=np.nanmax(abs(V))),"sky-mean LOS relative velocity (km/s, symmetric log)",0,f"heat_pm_mean_D_s{suf}.png",f"Collapse, point mass 7.5×10⁵² kg — {lab}\nsky-mean relative velocity (speeds scale as √M)")):
        f,ax=plt.subplots(figsize=(8.5,5.4)); ax.set_facecolor("0.85")
        pc=ax.pcolormesh(DD,S,Z,cmap="RdBu_r",norm=norm,shading="auto")
        if np.nanmin(Z)<lev<np.nanmax(Z): cs=ax.contour(DD,S,Z,levels=[lev],colors="k",linewidths=2.2); ax.clabel(cs,fmt=lambda v:"50%" if lev==50 else "0 km/s",fontsize=10)
        ax.axvline(58.2,color="k",lw=.8,ls=":"); ax.text(58.8,29,"paper D",fontsize=9)
        ax.set_xlabel("observer distance to CoE, D (Gpc)",fontsize=11); ax.set_ylabel("galaxy distance s (Gpc)",fontsize=11); ax.set_title(tt,fontsize=11.5)
        cb=f.colorbar(pc,ax=ax); cb.set_label(cl,fontsize=10.5); f.savefig(fl,dpi=130,bbox_inches="tight"); plt.close(f)
