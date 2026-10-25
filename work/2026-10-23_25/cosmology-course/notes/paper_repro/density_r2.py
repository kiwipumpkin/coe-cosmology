import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
c=1.0; D=lambda t:3.0+t
def bands(t,K):
    L=np.linspace(D(t)+c*t,D(t)-c*t,2*K*t+1); ro,ri=L[:-1],L[1:]
    return ro,ri,(ro+ri)/2,(4/3)*np.pi*(ro**3-ri**3)
res={}
for K in (1,2):
    ro,ri,m,V=bands(3,K)
    f=lambda p:(m[-1]**p/V[-1])/(m[0]**p/V[0])-1.03
    p=brentq(f,0.5,4); rho=m**p/V; rho0=rho[-1]
    assert abs(rho[-1]/rho[0]-1.03)<1e-10 and np.all(np.diff(rho)>0)   # rises inward = falls outward
    for t in (1,2,3):
        r=bands(t,K); rr=r[2]**p/r[3]; assert np.all(np.diff(rr)>0)
    res[K]=(p,rho0)
    print(f"{6*K} shells: p={p:.4f} rho(D-R)={rho[-1]:.6e} rho(D+R)={rho[0]:.6e} ratio={rho[-1]/rho[0]:.6f} content out/in={(m[0]/m[-1])**p:.4f}")
np.save("p_r2.npy",[res[1][0],res[1][1],res[2][0],res[2][1]])
fig,ax=plt.subplots(figsize=(11,6.5))
for K,ls in ((1,"-"),(2,"--")):
    p,rho0=res[K]
    for t,col in ((1,"tab:blue"),(2,"tab:green"),(3,"tab:red")):
        ro,ri,m,V=bands(t,K); rho=m**p/V/rho0
        ax.step(np.r_[ri[::-1],ro[0]],np.r_[rho[::-1],rho[0]],where="post",color=col,ls=ls,lw=2.2 if K==1 else 1.4,
                label=f"{6*K if t==3 else 2*K*t} shells, t={t} (p={p:.3f})")
for t in (1,2,3): ax.axvline(D(t)+t,color="gray",ls=":"); ax.text(D(t)+t,1.001,f" D+R (t={t})",fontsize=9)
ax.axvline(3,color="k",ls=":"); ax.text(3,1.001," D−R (all t)",fontsize=9)
ax.set_xlabel("distance from CoE (c·Δt)"); ax.set_ylabel("density ÷ density of innermost band at t=3")
ax.set_title("Content per shell ∝ r^p, p set so density at D−R is 3% above D+R at t=3\nsolid = shells c·Δt wide, dashed = shells c·Δt/2 wide")
ax.grid(alpha=.3); ax.legend(fontsize=9); plt.tight_layout(); plt.savefig("density_r2.png",dpi=200)
