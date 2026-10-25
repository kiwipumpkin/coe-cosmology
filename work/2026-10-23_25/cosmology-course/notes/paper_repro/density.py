import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
c=1.0; D=lambda t:3.0+t
def bands(t):
    L=np.linspace(D(t)+c*t,D(t)-c*t,2*t+1); ro,ri=L[:-1],L[1:]
    return ro,ri,(ro+ri)/2,(4/3)*np.pi*(ro**3-ri**3)   # index 0 = shell 1 (outermost)
ro,ri,m3,V3=bands(3)
# rho_in/rho_out = (V_out/V_in)*exp(-(r_in-r_out)/L) = 1.03
Lsc=(m3[0]-m3[-1])/np.log(1.03*V3[-1]/V3[0])
rref=m3[-1]; mass=lambda r: np.exp(-(r-rref)/Lsc)   # m0=1 at innermost band of t=3
rat=(mass(m3[-1])/V3[-1])/(mass(m3[0])/V3[0]); assert abs(rat-1.03)<1e-12
print("L=",Lsc,"ratio",rat,"content outer/inner",mass(m3[0]))
fig,ax=plt.subplots(1,2,figsize=(16,6.5),sharex=True)
cols={1:"tab:blue",2:"tab:green",3:"tab:red"}
for t in (1,2,3):
    ro,ri,m,V=bands(t)
    for a,rho,ttl in ((ax[0],1/V,"equal content per shell"),(ax[1],mass(m)/V,f"3% lock (content ∝ exp(+(r−{rref})/{-Lsc:.2f}), grows outward)")):
        a.step(np.r_[ri[::-1],ro[0]],np.r_[rho[::-1],rho[0]],where="post",color=cols[t],lw=2,label=f"t={t}")
        a.axvline(D(t)-t,color=cols[t],ls=":"); a.axvline(D(t)+t,color=cols[t],ls="--")
        a.set_title(ttl); a.set_xlabel("distance from CoE (c·Δt)"); a.set_ylabel("density (content ÷ shell volume)")
        if t==3: print(ttl,[f"{x:.5f}" for x in rho[::-1]])
for a in ax:
    a.text(3,a.get_ylim()[1]*0.95," D−R (all t)",fontsize=10)
    a.legend(title="dotted = D−R, dashed = D+R"); a.grid(alpha=.3)
plt.tight_layout(); plt.savefig("density_vs_distance.png",dpi=200)
np.save("L3pct.npy",[Lsc,rref])
