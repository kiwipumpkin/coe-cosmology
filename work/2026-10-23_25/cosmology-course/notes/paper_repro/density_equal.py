import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
c=1.0; D=lambda t:3.0+t
def bands(t,K):
    L=np.linspace(D(t)+c*t,D(t)-c*t,2*K*t+1); ro,ri=L[:-1],L[1:]
    return ro,ri,(ro+ri)/2,(4/3)*np.pi*(ro**3-ri**3)
fig,ax=plt.subplots(1,2,figsize=(18,7))
cols={1:"tab:blue",2:"tab:green",3:"tab:red"}
for K,ls,lw in ((1,"-",2.4),(2,"--",1.4)):
    ro,ri,m,V=bands(3,K); ref=1/V[-1]
    p=brentq(lambda p:(m[-1]**p/V[-1])/(m[0]**p/V[0])-1.03,0.5,4)
    for t in (1,2,3):
        ro_,ri_,mm,VV=bands(t,K); rho=(1/VV)/ref; assert np.all(np.diff(rho)>0)
        x=np.r_[ri_[::-1],ro_[0]]
        ax[0].step(x,np.r_[rho[::-1],rho[0]],where="post",color=cols[t],ls=ls,lw=lw,label=f"{2*K*t} shells, t={t}")
    print(f"{6*K} shells equal: rho(D-R)={1/V[-1]:.5e} rho(D+R)={1/V[0]:.5e} ratio={V[0]/V[-1]:.4f}")
    x=np.r_[ri[::-1],ro[0]]
    for q,lab,col in ((np.ones_like(m),"equal content",'tab:red'),(m**p,f"3% lock (r^{p:.3f})",'tab:purple'),(V,"flat (∝ volume)",'tab:gray')):
        rho=q/V; rho=rho/rho[-1]
        ax[1].step(x,np.r_[rho[::-1],rho[0]],where="post",color=col,ls=ls,lw=lw,label=f"{lab}, {6*K} shells")
for a in ax:
    a.axvline(3,color="k",ls=":"); a.axvline(9,color="k",ls=":")
    a.set_xlabel("distance from CoE (c·Δt)"); a.grid(alpha=.3,which="both")
for t in (1,2): ax[0].axvline(D(t)+t,color="gray",ls=":")
ax[0].set_title("Equal matter+light per shell: density = 1 ÷ shell volume\nsolid = cΔt wide, dashed = cΔt/2 wide")
ax[0].set_ylabel("density ÷ innermost band at t=3"); ax[0].legend(fontsize=9)
ax[1].set_yscale("log"); ax[1].set_title("t = 3: equal content vs 3% lock vs flat (log scale)")
ax[1].set_ylabel("density ÷ innermost band"); ax[1].legend(fontsize=9)
for a in ax:
    y=a.get_ylim()[1]; a.text(3,y," D−R",va="top"); a.text(9,y,"D+R (t=3) ",va="top",ha="right")
plt.tight_layout(); plt.savefig("density_equal.png",dpi=180)
