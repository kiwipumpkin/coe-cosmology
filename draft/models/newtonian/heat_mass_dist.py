"""Collapse, kappa=1.05, m=-4 (collapse_ls family). Point mass 7.5e52 kg fixed + uniform shell (43.9-72.5 Gpc) with density
f x 2.9e-27 kg/m^3. Pure point mass: mass only scales speeds (fractions fixed). Adding shell mass changes M(r) with radius, so the shape changes."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm,SymLogNorm
G=6.674e-11; GPC=3.0857e25; MP=7.52e52; RHO=2.9e-27; D=58.2; kap,m=1.05,-4
th=np.radians(np.linspace(0,180,721)); ct,st=np.cos(th),np.sin(th); w=st/st.sum()
a,b=43.9,72.5
def stats(f,s):
    def M(r): rr=np.clip(r,a,b); return MP+f*4/3*np.pi*RHO*GPC**3*(rr**3-a**3)
    def v(r): r=np.asarray(r,float); rmax=kap*D*(r/D)**m; return -np.sqrt(np.maximum(2*G*M(r)*(1/r-1/rmax)/GPC,0))
    x=D+s*ct; y=s*st; r=np.hypot(x,y); vl=v(r)*(x*ct+y*st)/r-v(np.array(D))*ct
    return (w*(vl>0)).sum(),(w*vl).sum()/1e3
F=np.logspace(-5,0,61); S=np.linspace(1,30,59)
R=np.zeros((len(S),len(F))); V=np.zeros_like(R)
for i,s in enumerate(S):
    for j,f in enumerate(F): R[i,j],V[i,j]=stats(f,s)
Mshell=4/3*np.pi*RHO*GPC**3*(b**3-a**3)
print("shell mass at f=1: %.3g kg; equals point mass at f=%.2g"%(Mshell,MP/Mshell))
print("frac range %.3f-%.3f; mean v range %.0f to %.0f km/s"%(R.min(),R.max(),V.min(),V.max()))
for s in (1,8,14.3,30):
    i=np.argmin(abs(S-s)); r=R[i]; v=V[i]
    print(f"s={S[i]:.1f}: frac f=1e-5 {r[0]:.2f}, f=1 {r[-1]:.2f}; mean v {v[0]:.0f} -> {v[-1]:.0f}; 50% cross at f={F[np.where(np.diff(np.sign(r-.5)))[0]] if np.any(np.diff(np.sign(r-.5))) else 'none'}; v=0 cross at f={F[np.where(np.diff(np.sign(v)))[0]] if np.any(np.diff(np.sign(v))) else 'none'}")
for Z,cmap_norm,lab,lev,fn,title in ((R*100,TwoSlopeNorm(50,0,100),"sky fraction receding (%)",[50],"heat_frac_mass_dist.png","Collapse (κ = 1.05, m = −4): receding sky fraction"),
                                    (V,SymLogNorm(linthresh=1000,vmin=-abs(V).max(),vmax=abs(V).max(),base=10),"sky-mean LOS relative velocity (km/s, symmetric log scale)",[0],"heat_mean_mass_dist.png","Collapse (κ = 1.05, m = −4): sky-mean relative velocity")):
    f,ax=plt.subplots(figsize=(8.5,5.4))
    pc=ax.pcolormesh(F,S,Z,cmap="RdBu_r",norm=cmap_norm,shading="auto"); ax.set_xscale("log")
    if Z.min()<lev[0]<Z.max(): cs=ax.contour(F,S,Z,levels=lev,colors="k",linewidths=2.2); ax.clabel(cs,fmt=lambda v:("50%" if "frac" in fn else "0 km/s"),fontsize=10)
    sec=ax.secondary_xaxis("top",functions=(lambda x:x*Mshell/MP,lambda y:y*MP/Mshell)); sec.set_xlabel("shell mass / central point mass",fontsize=10)
    ax.set_xlabel("shell density multiple  (× 2.9×10⁻²⁷ kg/m³; central point mass 7.5×10⁵² kg fixed)",fontsize=10.5)
    ax.set_ylabel("galaxy distance s (Gpc)",fontsize=11); ax.set_title(title,fontsize=12,pad=28)
    cb=f.colorbar(pc,ax=ax); cb.set_label(lab,fontsize=10.5); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
