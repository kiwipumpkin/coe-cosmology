import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
c=1.0; D=lambda t:3.0+t
def bands(t,K):
    L=np.linspace(D(t)+c*t,D(t)-c*t,2*K*t+1); ro,ri=L[:-1],L[1:]
    assert abs(ro[0]-(D(t)+c*t))<1e-9 and abs(ri[-1]-(D(t)-c*t))<1e-9
    return ro,ri,(ro+ri)/2,(4/3)*np.pi*(ro**3-ri**3)   # idx0 = shell 1 (outermost)
fig,ax=plt.subplots(2,2,figsize=(16,11),sharex=True)
cols={1:"tab:blue",2:"tab:green",3:"tab:red"}
for K,ls,lw in ((1,"-",2.4),(2,"--",1.4)):
    ro,ri,m,V=bands(3,K)
    p=brentq(lambda p:(m[-1]**p/V[-1])/(m[0]**p/V[0])-1.03,0.5,4)
    laws={"A":lambda m,V: m**p, "B":lambda m,V: V}
    refs={k:laws[k](m,V)[-1] for k in laws}             # innermost shell at t=3
    for t in (1,2,3):
        ro,ri,mm,VV=bands(t,K)
        for j,k in enumerate("AB"):
            q=laws[k](mm,VV); rho=q/VV
            if k=="A": assert np.all(np.diff(rho)>0)
            else: assert np.ptp(rho)<1e-12*rho.max()
            if t==3 and k=="A": assert abs(rho[-1]/rho[0]-1.03)<1e-10
            rn=rho/(refs[k]/V[-1]); qn=q/refs[k]
            x=np.r_[ri[::-1],ro[0]]
            ax[0,j].step(x,np.r_[rn[::-1],rn[0]],where="post",color=cols[t],ls=ls,lw=lw,label=f"{2*K*t} shells, t={t}")
            ax[1,j].step(x,np.r_[qn[::-1],qn[0]],where="post",color=cols[t],ls=ls,lw=lw,label=f"{2*K*t} shells, t={t}")
            if t==3: print(f"{6*K} shells {k}: rho(D-R)={rho[-1]:.5e} rho(D+R)={rho[0]:.5e} ratio={rho[-1]/rho[0]:.6f} content out/in={q[0]/q[-1]:.4f}"+(f" p={p:.4f}" if k=="A" else ""))
ax[0,0].set_title("A: 3% lock — content per shell ∝ r_mid^p (p: 6 shells 1.960, 12 shells 1.968)")
ax[0,1].set_title("B: flat density — content per shell ∝ shell volume")
ax[0,0].set_ylabel("density ÷ innermost band at t=3"); ax[1,0].set_ylabel("content (matter+light) per shell ÷ innermost shell at t=3")
ax[0,1].set_ylim(0.97,1.01)
for a in ax.flat:
    a.axvline(3,color="k",ls=":"); a.axvline(9,color="k",ls=":")
    a.text(3,a.get_ylim()[1]," D−R",va="top"); a.text(9,a.get_ylim()[1]," D+R (t=3)",va="top",ha="right")
    a.grid(alpha=.3); a.legend(fontsize=8,title="solid = cΔt wide, dashed = cΔt/2 wide",title_fontsize=8)
for a in ax[1]: a.set_xlabel("distance from CoE (c·Δt)")
plt.tight_layout(); plt.savefig("density_and_mass.png",dpi=180)
