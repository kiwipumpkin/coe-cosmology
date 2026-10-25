import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import quad
a,b=3.0,9.0; k=np.log(1.03)/np.log(3)          # power-law exponent 0.0269
Vtot=4/3*np.pi*(b**3-a**3); rbar=1/Vtot          # mean density in window (mass 1)
def norm(f):
    Z=quad(lambda r:4*np.pi*r*r*f(r),a,b)[0]; return lambda r: f(r)/Z
fam=[]
fam.append(("power law ρ ∝ r^−%.4f"%k,"tab:purple",norm(lambda r:r**-k)))
for n in (1,2,3):
    BA=0.03/(3.0**-n-1.03*9.0**-n); fam.append((f"A + B/r^{n}  (B/A={BA:.4f})",plt.cm.Reds(0.4+0.2*n),norm(lambda r,n=n,BA=BA:1+BA/r**n)))
for L in (1.0,3.0,10.0):
    BA=0.03/(np.exp(-3/L)-1.03*np.exp(-9/L)); fam.append((f"A + B·exp(−r/{L:g})  (B/A={BA:.4f})",plt.cm.Greens(0.35+0.05*L),norm(lambda r,L=L,BA=BA:1+BA*np.exp(-r/L))))
for cA in (-6e-4,-3e-4,0.0,3e-4,6e-4):
    bA=-0.03*(1+9*cA)/(3*2.03)
    fam.append((f"A[1 + β(r−6) + γ(r−6)²]  (γ={cA:+.0e}, β={bA:.5f})",plt.cm.Blues(0.4+600*abs(cA)+(0.1 if cA>0 else 0)),norm(lambda r,cA=cA,bA=bA:1+bA*(r-6)+cA*(r-6)**2)))
pl=fam[0][2]; rw=np.linspace(a,b,601)
rows=[]
for name,col,f in fam:
    assert abs(f(3)/f(9)-1.03)<1e-12 and abs(quad(lambda r:4*np.pi*r*r*f(r),a,b)[0]-1)<1e-10
    dev=np.max(np.abs(f(rw)/pl(rw)-1)); rows.append((name,dev,[f(x)/rbar for x in (0.5,3,9,12)]))
    assert dev<1.5e-2, (name,dev)
fig,ax=plt.subplots(3,1,figsize=(15,16))
rin=np.linspace(a,b,300); rlo=np.linspace(0.3,a,300); rhi=np.linspace(b,12,200)
for name,col,f in fam:
    lw=3 if name.startswith("power") else 1.6
    for j in (0,1):
        ax[j].plot(rin,f(rin)/rbar,color=col,lw=lw,label=name if j==0 else None)
    for rr in (rlo,rhi): ax[0].plot(rr,f(rr)/rbar,color=col,lw=lw*0.8,ls="--")
    g=lambda r:4*np.pi*r*r*f(r)
    ax[2].plot(rin,g(rin),color=col,lw=lw)
    for rr in (rlo,rhi): ax[2].plot(rr,g(rr),color=col,lw=lw*0.8,ls="--")
ax[0].set_ylim(0.9,1.4); ax[1].set_xlim(a,b)
ax[0].set_title("density ÷ mean density in 3–9 (y clipped at 1.4; steep inner laws rise off-scale toward CoE)")
ax[1].set_title("zoom 3–9: curves nearly overlap (max deviation from power law 0.09%–0.63%, see report)")
ax[2].set_title("mass per unit width, 4πr²ρ (total = 1 over 3–9)")
for j,a_ in enumerate(ax):
    a_.grid(alpha=.3); a_.set_xlabel("distance from CoE (c·Δt)")
    if j!=1:
        a_.set_xlim(0,12); a_.axvspan(a,b,color="khaki",alpha=0.25,zorder=-1)
        a_.set_xticks([0,3,6,9,12]); a_.set_xticklabels(["CoE\n0","D−R\n3","D\n6","D+R\n9","D+2R\n12"])
    else:
        a_.set_xticks([3,4,5,6,7,8,9]); a_.set_xticklabels(["D−R\n3","4","5","D\n6","7","8","D+R\n9"])
    for x in (3,6,9,12): a_.axvline(x,color="k",lw=0.8,ls=":")
ax[0].set_ylabel("ρ / mean ρ"); ax[1].set_ylabel("ρ / mean ρ"); ax[2].set_ylabel("mass per unit width")
fig.legend(*ax[0].get_legend_handles_labels(),loc="center left",bbox_to_anchor=(0.70,0.5),fontsize=9,frameon=False)
fig.suptitle("Density laws that all give ρ(3)/ρ(9) = 1.03 and total mass 1 in 3–9, differing outside. Solid in 3–9, dashed outside.",fontsize=12)
plt.tight_layout(rect=[0,0,0.70,0.96]); plt.savefig("level_sets_3pct.png",dpi=160)
for nm,d,v in rows: print(f"{nm:55s} maxdev={d:.2e}  ρ/mean at 0.5,3,9,12 = "+", ".join(f"{x:.4f}" for x in v))
