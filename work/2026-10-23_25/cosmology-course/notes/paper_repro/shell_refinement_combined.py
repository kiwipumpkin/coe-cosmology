import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("shell_refinement.py").read().split("names=")[0])
names={"eq":"equal mass per shell","flat":"flat density","3pct":"3% density drop"}
col={"eq":"tab:red","flat":"tab:gray","3pct":"tab:purple"}
sty={6:(":",0.45),12:("--",0.65),24:("-",0.85)}
fig,AX=plt.subplots(2,2,figsize=(17,12)); ax=AX[0]
for k in names:
    g=cont(k)
    for N in Ns:
        e,m,mw,rho,p=stairs(k,N); ls,al=sty[N]
        lab=f"{names[k]}, {N} shells"+(f" (p={p:.4f})" if p else "")
        W=(3.2 if k=='flat' else 1.4); Z=(1 if k=='flat' else 3)
        ax[0].stairs(mw,e,color=col[k],ls=ls,alpha=al,lw=W,zorder=Z,baseline=None,label=lab)
        ax[1].stairs(rho,e,color=col[k],ls=ls,alpha=al,lw=W,zorder=Z,baseline=None)
    ax[0].plot(r,g(r),color=col[k],lw=6 if k=='flat' else 2.5,zorder=2 if k=='flat' else 4,label=f"{names[k]}, continuum"+(f" (p={qc:.4f})" if k=="3pct" else ""))
    ax[1].plot(r,g(r)/(4*np.pi*r**2),color=col[k],lw=6 if k=='flat' else 2.5,zorder=2 if k=='flat' else 4)
gF,g3=cont("flat"),cont("3pct")
for N in Ns:
    e,_,mwF,rF,_=stairs("flat",N); e,_,mw3,r3,_=stairs("3pct",N); ls,al=sty[N]
    AX[1,0].stairs(mw3/mwF,e,color="tab:purple",ls=ls,alpha=al,lw=1.6,baseline=None)
    AX[1,1].stairs(r3/rF,e,color="tab:purple",ls=ls,alpha=al,lw=1.6,baseline=None)
AX[1,0].plot(r,g3(r)/gF(r),color="tab:purple",lw=2.5); AX[1,1].plot(r,g3(r)/gF(r),color="tab:purple",lw=2.5)
for a_ in AX[1]: a_.axhline(1,color="tab:gray",lw=4,zorder=0)
AX[1,1].set_ylim(0.96,1.03)
AX[1,0].set_title("zoom: mass per width, 3% ÷ flat"); AX[1,1].set_title("zoom: density, 3% ÷ flat (gray line = flat)")
AX[1,0].set_ylabel("ratio to flat"); AX[1,1].set_ylabel("ratio to flat")
ax[0].set_title("mass per unit width (mass per shell ÷ shell width)"); ax[1].set_title("density (mass ÷ shell volume)")
ax[0].set_ylabel("mass per unit width"); ax[1].set_ylabel("density")
for a_ in AX.flat:
    a_.set_xlim(3,9); a_.grid(alpha=.3); a_.set_xlabel("distance from CoE (c·Δt)")
    a_.set_xticks(range(3,10)); a_.set_xticklabels(["D−R\n3","4","5","6","7","8","D+R\n9"])
fig.legend(*ax[0].get_legend_handles_labels(),loc="center left",bbox_to_anchor=(0.80,0.5),fontsize=9,frameon=False)
fig.suptitle("t = 3, exactly-c geometry. Total mass between D−R and D+R = 1 in every case. Dotted 6, dashed 12, thin solid 24 shells; thick = continuum",fontsize=12)
plt.tight_layout(rect=[0,0,0.80,0.94]); plt.savefig("shell_refinement_combined.png",dpi=170)
