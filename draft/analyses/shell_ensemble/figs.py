import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":16,"axes.titlesize":17,"axes.labelsize":17,"xtick.labelsize":15,"ytick.labelsize":15,"legend.fontsize":14})
o=json.load(open("ensemble.json")); F=1e-5
LAB={"a":r"(a) $10^{48}$–$10^{49}$ kg","b":r"(b) $10^{49.5}$–$10^{50.5}$ kg","c":r"(c) $10^{51}$–$10^{52}$ kg"}
COL={"a":"C2","b":"C1","c":"C3"}
for tag,title in (("bare","objects anywhere beyond 14.3 Gpc"),("buf","whole object beyond 14.3 Gpc")):
    fig,ax=plt.subplots(2,2,figsize=(18,13),sharex=True)
    bins=np.logspace(np.log10(0.8),np.log10(5e3),140)
    for i,N in enumerate((200,500)):
        for j,(key,lbl) in enumerate((("rms","random web as drawn"),("rms_best","best global amplitude per realization"))):
            a=ax[i,j]
            for lev in "abc":
                d=np.load(f"raw_{tag}_N{N}_{lev}.npz")[key]/F
                a.hist(d,bins=bins,histtype="stepfilled",alpha=0.45,color=COL[lev],label=LAB[lev]+f", median {np.median(d):.0f}×")
            for x,ls,c,t in ((237.26,"-","k","no web 237×"),(83.4,"--","C0","four clumps 83×"),(1,":","purple","floor 1×")):
                a.axvline(x,ls=ls,color=c,lw=2.5,label=t)
            a.set_xscale("log"); a.set_yscale("log"); a.set_title(f"N = {N}: {lbl}"); a.grid(alpha=.3)
            a.set_ylabel("realizations (of 1000)"); a.legend(loc="upper left",fontsize=12.5,framealpha=.9)
            if i==1: a.set_xlabel("residual RMS ÷ $10^{-5}$ floor")
    fig.suptitle(f"Random cluster+void web in the CoE shell (43.9–72.5 Gpc), zero net mass — {title}",fontsize=19)
    fig.tight_layout(rect=(0,0,1,.96)); fig.savefig(f"hist_rms_{tag}.png",dpi=110); plt.close(fig)

# example realization
ex=o["examples"]["bare_N500_c"]; P=np.array(ex["P"]); M=np.array(ex["M"]); ell=np.radians(o["ell_coe"])
C=58.2*np.array([np.cos(ell),np.sin(ell)])
sz=40*np.abs(M)/np.abs(M).max()*8+6
fig=plt.figure(figsize=(21,11))
a=fig.add_subplot(1,2,1); t=np.linspace(0,2*np.pi,400)
for r,ls in ((43.9,"--"),(72.5,"--")): a.plot(C[0]+r*np.cos(t),C[1]+r*np.sin(t),"k",ls=ls,lw=1.5)
a.fill(14.3*np.cos(t),14.3*np.sin(t),color="C0",alpha=.25,label="vision sphere 14.3 Gpc")
cl=M>0
a.scatter(P[cl,0],P[cl,1],s=sz[cl],c="red",alpha=.6,edgecolor="k",lw=.3,label="clusters (δ>0)")
a.scatter(P[~cl,0],P[~cl,1],s=sz[~cl],c="blue",alpha=.6,edgecolor="k",lw=.3,label="voids (δ<0)")
a.plot([-30*np.cos(ell),C[0]*1.0],[-30*np.sin(ell),C[1]],"k-",lw=1.2); a.scatter(*C,marker="x",s=250,c="k",lw=3,label="CoE (58.2 Gpc)")
a.scatter(0,0,marker="*",s=400,c="gold",edgecolor="k",zorder=5,label="observer")
a.annotate(f"ℓ = {np.degrees(ell)%360:.1f}° (toward CoE)",xy=C*0.45,xytext=(-55,8),fontsize=14)
a.set_aspect("equal"); a.set_xlabel("x (Gpc), toward ℓ = 0°"); a.set_ylabel("y (Gpc), toward ℓ = 90°")
a.set_title("Top-down projection on the ring plane\n(dashed: shell radii in-plane; off-plane objects project inside)",fontsize=15); a.legend(loc="upper left",fontsize=12.5,framealpha=.9); a.grid(alpha=.3)
b=fig.add_subplot(1,2,2,projection="3d")
b.scatter(P[cl,0],P[cl,1],P[cl,2],s=sz[cl]*.6,c="red",alpha=.5)
b.scatter(P[~cl,0],P[~cl,1],P[~cl,2],s=sz[~cl]*.6,c="blue",alpha=.5)
u,v=np.mgrid[0:2*np.pi:30j,0:np.pi:15j]
b.plot_wireframe(14.3*np.cos(u)*np.sin(v),14.3*np.sin(u)*np.sin(v),14.3*np.cos(v),color="C0",lw=.5,alpha=.6)
b.plot_wireframe(C[0]+72.5*np.cos(u)*np.sin(v),C[1]+72.5*np.sin(u)*np.sin(v),72.5*np.cos(v),color="0.6",lw=.3,alpha=.25)
b.scatter(0,0,0,marker="*",s=300,c="gold",edgecolor="k"); b.scatter(C[0],C[1],0,marker="x",s=200,c="k")
b.plot([0,C[0]],[0,C[1]],[0,0],"k-",lw=1.5)
b.set_xlabel("x (Gpc)",labelpad=12); b.set_ylabel("y (Gpc)",labelpad=12); b.set_zlabel("z (Gpc)",labelpad=12)
b.set_title("3D view (blue mesh: vision sphere; grey: outer shell 72.5 Gpc)")
fig.subplots_adjust(); fig.suptitle(f"Example realization: N = 500, level (c) 10$^{{51}}$–10$^{{52}}$ kg, zero net mass (marker area ∝ |M|)",fontsize=19)
fig.tight_layout(rect=(0,0,1,.94)); fig.savefig("example_realization.png",dpi=110); plt.close(fig)
print("ok")
