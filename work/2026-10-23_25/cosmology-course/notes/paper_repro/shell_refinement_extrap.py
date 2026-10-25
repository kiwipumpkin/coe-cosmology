import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("shell_refinement.py").read().split("names=")[0])
names={"eq":"equal mass per shell","flat":"flat density","3pct":"3% density drop"}
col={"eq":"tab:red","flat":"tab:gray","3pct":"tab:purple"}; sty={6:(":",0.45),12:("--",0.65),24:("-",0.85)}
rin=np.linspace(3,9,400); rlo=np.linspace(0.02,3,300); rhi=np.linspace(9,12,200)
fig,AX=plt.subplots(2,2,figsize=(18,12)); ax=AX[0]
dens=lambda g: (lambda x: g(x)/(4*np.pi*x**2))
for k in names:
    g=cont(k); W=3.2 if k=="flat" else 1.4; Z=1 if k=="flat" else 3; LW=6 if k=="flat" else 2.5
    for N in Ns:
        e,m,mw,rho,p=stairs(k,N); ls,al=sty[N]
        ax[0].stairs(mw,e,color=col[k],ls=ls,alpha=al,lw=W,zorder=Z,baseline=None,label=f"{names[k]}, {N} shells"+(f" (p={p:.4f})" if p else ""))
        ax[1].stairs(rho,e,color=col[k],ls=ls,alpha=al,lw=W,zorder=Z,baseline=None)
    for j,f in enumerate((g,dens(g))):
        ax[j].plot(rin,f(rin),color=col[k],lw=LW,zorder=Z+1,label=(f"{names[k]}, continuum"+(f" p={qc:.4f}" if k=="3pct" else "")) if j==0 else None)
        for rr in (rlo,rhi): ax[j].plot(rr,f(rr),color=col[k],lw=LW*0.6,ls="--",zorder=Z+1)
    print(k,"mass/width at 0+,3,9,12:",[f"{g(x):.4g}" for x in (0.02,3,9,12)],"density:",[f"{dens(g)(x):.4g}" for x in (0.02,3,9,12)])
gF,g3=cont("flat"),cont("3pct")
for N in Ns:
    e,_,mwF,rF,_=stairs("flat",N); e,_,mw3,r3,_=stairs("3pct",N); ls,al=sty[N]
    for a_,v in ((AX[1,0],mw3/mwF),(AX[1,1],r3/rF)): a_.stairs(v,e,color="tab:purple",ls=ls,alpha=al,lw=1.6,baseline=None)
for a_ in AX[1]:
    a_.plot(rin,g3(rin)/gF(rin),color="tab:purple",lw=2.5)
    for rr in (rlo,rhi): a_.plot(rr,g3(rr)/gF(rr),color="tab:purple",lw=1.6,ls="--")
    a_.axhline(1,color="tab:gray",lw=4,zorder=0); a_.set_ylim(0.97,1.18)
print("3%/flat ratio at 0+,3,9,12:",[f"{g3(x)/gF(x):.4f}" for x in (0.02,3,9,12)])
ax[1].set_ylim(0,0.0025); ax[0].set_ylim(0,0.65)
ax[0].set_title("mass per unit width"); ax[1].set_title("density (y clipped at 0.0025; equal-mass ∝ 1/r² diverges at r→0)")
AX[1,0].set_title("zoom: mass per width, 3% ÷ flat"); AX[1,1].set_title("zoom: density, 3% ÷ flat (gray = flat)")
ax[0].set_ylabel("mass per unit width"); ax[1].set_ylabel("density"); AX[1,0].set_ylabel("ratio to flat"); AX[1,1].set_ylabel("ratio to flat")
for a_ in AX.flat:
    a_.set_xlim(0,12); a_.axvspan(3,9,color="khaki",alpha=0.25,zorder=-1); a_.grid(alpha=.3)
    for x in (3,6,9,12): a_.axvline(x,color="k",lw=1,ls=":")
    a_.set_xticks([0,3,6,9,12]); a_.set_xticklabels(["CoE\n0","D−R = R_rec\n3","D\n6","D+R\n9","D+2R\n12"]); a_.set_xlabel("distance from CoE (c·Δt)")
fig.legend(*ax[0].get_legend_handles_labels(),loc="center left",bbox_to_anchor=(0.80,0.5),fontsize=9,frameon=False)
fig.suptitle("t = 3, exactly-c geometry. Shaded = visible region 3–9 (total mass there = 1). Continuum solid in 3–9, dashed = extrapolated to 0 and to D+2R = 12.",fontsize=12)
plt.tight_layout(rect=[0,0,0.80,0.95]); plt.savefig("shell_refinement_extrap.png",dpi=170)
