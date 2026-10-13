"""Four parameters: D, R, M, v0. All shells share our energy: v(r)^2 = v0^2 + 2GM(1/r-1/D), sign of v0.
In units of v_esc(D)=sqrt(2GM/D): (v/v_esc)^2 = x^2 + D/r - 1, x=|v0|/v_esc, r/D = |(1,0)+(R/D) n|.
So LOS velocity / |v0| depends only on x and R/D (M, D, v0 enter only through x). Cell invalid (grey) if v^2<0 on any sightline."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
th=np.radians(np.linspace(0,180,721)); ct,st=np.cos(th),np.sin(th); w=st/st.sum()
def calc(x,q,sg):
    X=1+q*ct; Y=q*st; rho=np.hypot(X,Y); u2=x**2+1/rho-1
    if (u2<0).any() or x==0: return np.nan,np.nan
    vl=sg*(np.sqrt(u2)*(X*ct+Y*st)/rho - x*ct)/x       # in units of |v0|
    return (w*(vl>0)).sum(),(w*vl).sum()
# numeric check of the scaling claim
G=6.674e-11;GPC=3.0857e25
def phys(D,R,M,v0):
    X=D+R*ct;Y=R*st;r=np.hypot(X,Y); v=np.sign(v0)*np.sqrt(v0**2+2*G*M*(1/(r*GPC)-1/(D*GPC))); vl=v*(X*ct+Y*st)/r-v0*ct; return (w*(vl>0)).sum(),(w*vl).sum()/abs(v0)
a=phys(58.2,14.3,7.5e52,-1e7); vesc=np.sqrt(2*G*7.5e52/(58.2*GPC)); x=1e7/vesc
b=phys(20,20*14.3/58.2,7.5e52*20/58.2*4,-2e7)   # same x and R/D, different D,M,v0
print("check same x,R/D:",a,b,calc(x,14.3/58.2,-1))
Xs=np.linspace(0.02,3,150); Q=np.linspace(0.02,0.5,97)
out={}
for sg,nm in ((-1,"col"),(1,"exp")):
    F=np.full((len(Q),len(Xs)),np.nan); V=F.copy()
    for i,q in enumerate(Q):
        for j,x in enumerate(Xs): F[i,j],V[i,j]=calc(x,q,sg)
    out[nm]=(F,V)
    print(nm,"frac %.3f-%.3f mean %.3f..%.3f |v0|"%(np.nanmin(F),np.nanmax(F),np.nanmin(V),np.nanmax(V)))
    for q in (0.05,0.25,0.5):
        i=np.argmin(abs(Q-q)); f=F[i]; ok=~np.isnan(f)
        c=Xs[ok][np.where(np.diff(np.sign(f[ok]-.5)))[0]]; c0=Xs[ok][np.where(np.diff(np.sign(V[i][ok])))[0]]
        print(f"   R/D={Q[i]:.2f}: valid from x={Xs[ok][0]:.2f}; frac at x_min {f[ok][0]:.2f}, x=1 {f[np.argmin(abs(Xs-1))]:.2f}, x=3 {f[-1]:.2f}; 50% at x={np.round(c,2)}; mean0 at x={np.round(c0,2)}")
def hm(Z,norm,lev,lab,title,fn,cbl):
    f,ax=plt.subplots(figsize=(8.5,5.2)); ax.set_facecolor("0.85")
    pc=ax.pcolormesh(Xs,Q,Z,cmap="RdBu_r",norm=norm,shading="auto")
    if np.nanmin(Z)<lev<np.nanmax(Z): cs=ax.contour(Xs,Q,Z,levels=[lev],colors="k",linewidths=2.2); ax.clabel(cs,fmt=lambda v:lab,fontsize=10)
    ax.axvline(1,color="k",lw=.8,ls=":"); ax.text(1.02,0.48,"escape speed",fontsize=9)
    ax.set_xlabel("x = |v0| / v_esc(D)   (our speed / escape speed at our radius)",fontsize=10.5); ax.set_ylabel("R / D  (galaxy distance / our distance to CoE)",fontsize=10.5)
    ax.set_title(title,fontsize=11.5); cb=f.colorbar(pc,ax=ax); cb.set_label(cbl,fontsize=10.5)
    f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
for nm,word in (("col","Collapsing (v0 < 0)"),("exp","Expanding (v0 > 0)")):
    F,V=out[nm]
    hm(F*100,TwoSlopeNorm(50,0,100),50,"50%",f"{word}: receding sky fraction\n(grey = some shells bound beyond turnaround, v² < 0)",f"simple4_frac_{nm}.png","sky fraction receding (%)")
    m=np.nanmax(abs(V)); hm(V,TwoSlopeNorm(0,-m,m),0,"0",f"{word}: sky-mean LOS relative velocity / |v0|",f"simple4_mean_{nm}.png","sky-mean relative velocity  (units of |v0|)")
print("\nD=58.2 Gpc, M=7.5e52 kg: v_esc = %.0f km/s = %.4f c"%(vesc/1e3,vesc/2.998e8))
for x in (0.25,0.5,1,1.5,2,3): print(f"  x={x}: |v0| = {x*vesc/1e3:,.0f} km/s ({x*vesc/2.998e8:.3f} c)")
