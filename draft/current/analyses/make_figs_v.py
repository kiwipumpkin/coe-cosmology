import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
exec(open("scan2.py").read().split("out={")[0])  # R, th, phi, field, models
S=json.load(open("scan.json")); KG=7.7e51/0.0039; rows=[]
def clumps(fam,N,p):
    if fam=="sym":
        npair=N//2; L=[]; i=1
        for _ in range(npair): L+=[(phi+p[i+1],p[i]),(phi-p[i+1],p[i])]; i+=2
        for s in S["sym"][str(N)]["singles"]: L.append((s,p[i])); i+=1
        return L, sym_model(np.array(p),npair,S["sym"][str(N)]["singles"])
    return [(p[2+2*j],p[1+2*j]) for j in range(N)], gen_model(np.array(p),N)
plt.rcParams.update({"font.size":15}); paths=[]
for fam,name in [("sym","symmetric pairs about axis"),("gen","free angle per clump")]:
  for N in range(4,9):
    r=S[fam][str(N)]; p=r["p"]; D=p[0]; L,m=clumps(fam,N,p); left=R-m; c=rms(left)
    masses=[k*KG for _,k in L]; fv=100*(1-np.mean(left**2)/np.mean(R**2))
    phys=all(x>0 for x in masses) and max(abs(x) for x in masses)<1e53 and D<19.9
    rows.append((N,fam,c,D,D*14.3,masses,phys))
    fig=plt.figure(figsize=(18,21)); gs=GridSpec(3,2,figure=fig,hspace=0.72,wspace=0.30,top=0.92,bottom=0.06,left=0.07,right=0.98)
    a=fig.add_subplot(gs[0,0]); a.plot(ELL,R*1e3,"C3",lw=2.5,label="residual after special relativity (SR) + mass")
    a.plot(ELL,m*1e3,"k--",lw=2,label=f"{N} clumps, tidal pull"); a.set_xlim(0,360); a.grid(alpha=.3)
    a.set_xlabel(r"Galactic longitude $\ell$ (deg)"); a.set_ylabel(r"$\Delta T/T\times10^{3}$")
    a.legend(loc="lower center",bbox_to_anchor=(0.5,1.01),frameon=False,fontsize=14)
    a.set_title(f"{N} clumps ({name}) explain {fv:.2f}% of post-SR variance",pad=72,fontsize=14)
    b=fig.add_subplot(gs[0,1],projection="polar"); rmax=max(2.35,D*1.25)
    b.plot(np.linspace(0,2*np.pi,400),np.ones(400),"C0",lw=2,label="visible universe, 14.3 Gpc")
    for (A,k),M in zip(L,masses):
        b.scatter(np.radians(A),D,s=400,color="C3" if M>0 else "C2",edgecolor="k",zorder=3)
        b.text(np.radians(A),D+0.13*rmax,f"{A%360:.0f}°\n{M:.1e} kg",ha="center",va="center",fontsize=11)
    b.scatter([],[],s=200,color="C3",label="mass concentration"); b.scatter([],[],s=200,color="C2",label="void")
    b.plot([np.radians(phi)]*2,[0,rmax*.95],"k",lw=1.5); b.plot([np.radians(phi+180)]*2,[0,rmax*.95],"k--",lw=1.5)
    b.set_ylim(0,rmax); b.set_yticklabels([]); b.legend(loc="lower center",bbox_to_anchor=(0.5,1.06),frameon=False,fontsize=13,ncol=2)
    b.set_title(f"clumps at {D*14.3:.1f} Gpc"+("  (D at fit bound!)" if D>19.9 else ""),pad=95,fontsize=15)
    e=fig.add_subplot(gs[1,:]); e.plot(ELL,left*1e3,"C0",lw=2.4,label=f"leftover after {N} clumps (residual − model)")
    e.axhline(0,color="k",lw=.7); e.axhline(1e-2,ls="--",color="purple",lw=1.8); e.axhline(-1e-2,ls="--",color="purple",lw=1.8,label=r"CMB level $\pm10^{-5}$")
    e.set_xlim(0,360); e.grid(alpha=.3); e.set_xlabel(r"Galactic longitude $\ell$ (deg)"); e.set_ylabel(r"$\Delta T/T\times10^{3}$")
    e.legend(loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=2,frameon=False,fontsize=14)
    e.set_title(f"Residual after {N} clumps: RMS = {c:.2e} (max |ΔT/T| = {np.abs(left).max():.2e}); top modes: {', '.join(modes(left))}",pad=52,fontsize=14)
    g_=fig.add_subplot(gs[2,:]); ms=np.arange(1,16)
    amp=lambda x:[np.hypot(2*np.mean(x*np.cos(k*th)),2*np.mean(x*np.sin(k*th))) for k in ms]
    g_.bar(ms-0.2,np.array(amp(R))*1e3,0.4,color="C3",label="post-SR residual"); g_.bar(ms+0.2,np.array(amp(left))*1e3,0.4,color="C0",label=f"leftover after {N} clumps")
    g_.set_yscale("log"); g_.set_xticks(ms); g_.set_xlabel("multipole m"); g_.set_ylabel(r"amplitude $\times10^{3}$"); g_.axhline(1e-2,ls="--",color="purple",label=r"$10^{-5}$")
    g_.grid(alpha=.3,axis="y"); g_.legend(loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=3,frameon=False,fontsize=14)
    f=f"figs_v/clumps_{fam}_N{N}.png"; fig.savefig(f,dpi=80,bbox_inches="tight"); plt.close(fig); paths.append(f)
pos=[f for f,(N,fam,c,D,G,M,ph) in zip(paths,rows) if all(x>0 for x in M)]
neg=[f for f,(N,fam,c,D,G,M,ph) in zip(paths,rows) if not all(x>0 for x in M)]
for name,L in [("clumps_sheet_masses_only.png",pos),("clumps_sheet_with_voids.png",neg)]:
    n=len(L); cols=2 if n>1 else 1; rws=(n+cols-1)//cols
    fig,ax=plt.subplots(rws,cols,figsize=(12*cols,13.2*rws),squeeze=False)
    for A in ax.flat: A.axis("off")
    for i,f in enumerate(L): ax.flat[i].imshow(plt.imread(f))
    fig.tight_layout(); fig.savefig(name,dpi=60); plt.close(fig); print(name,L)
