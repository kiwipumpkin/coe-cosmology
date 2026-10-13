"""Generate one network, measure the rule statistics, compare with targets -> stats_check.png, stats_check.json"""
import os,sys,json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ME=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,ME); import vmn_gen as G, lss
Mmin=float(sys.argv[1]) if len(sys.argv)>1 else 3e15
g=G.generate(Mmin=Mmin,seed=7); p=g["params"]; vi=g["void_in"]; L=g["links"]; R=p["R_out"]
# per-link geometry
dvec=g["XC"][L[:,0]]-g["S"][L[:,1]]; dist=np.linalg.norm(dvec,axis=1); dRv=dist/g["Reff"][L[:,1]]
# angular spread: resultant length of unit vectors to a void's clusters, vs isotropic with the same k
uv=dvec/dist[:,None]; Rsum=np.zeros((len(g["S"]),3))
for k in range(3): Rsum[:,k]=np.bincount(L[:,1],weights=uv[:,k],minlength=len(g["S"]))
kv=g["kv"]; sel=vi&(kv>=2); Rbar=np.linalg.norm(Rsum[sel],axis=1)/kv[sel]
rng=np.random.default_rng(3); Riso=[]
for k in kv[sel][:20000]:
    u=rng.normal(size=(k,3)); u/=np.linalg.norm(u,axis=1)[:,None]; Riso.append(np.linalg.norm(u.sum(0))/k)
Riso=np.array(Riso)
# mass vs mean bordering void radius
mR=np.bincount(L[:,0],weights=g["Reff"][L[:,1]],minlength=len(g["MC"]))/g["kc"]
slope=np.polyfit(np.log10(mR),np.log10(g["MC"]),1)[0]
# voids per cluster vs mass
mb=np.logspace(np.log10(Mmin),np.log10(g["MC"].max()),8); kcm=[g["kc"][(g["MC"]>=a)&(g["MC"]<b)].mean() for a,b in zip(mb[:-1],mb[1:])]
# compensation and budget (inside R)
V=4/3*np.pi*R**3; fvol=(p["fV"]*g["Vcell"][vi]).sum()/V
net=(-g["D"]+g["W"]+g["share"])[vi]/(lss.RHO_M*g["Vcell"][vi]); sig=lss.sigma(g["Reff"][vi])
clmass=g["MC"].sum()/(lss.RHO_M*V)
lm=np.linspace(np.log(Mmin),np.log(1e18),3000); clmass_t=np.trapezoid(lss.cluster_dn_dlnM(np.exp(lm))*np.exp(lm),lm)/lss.RHO_M
st=dict(Mmin=Mmin,n_voids_inside=int(vi.sum()),n_clusters=int(len(g["MC"])),merge_lam=g["lam"],spacing_Mpc=g["spacing"],
  clusters_per_void_mean=float(kv[vi].mean()),clusters_per_void_target=p["kv_mean"],
  voids_per_cluster_mean=float(g["kc"].mean()),voids_per_cluster_target=p["kc_target"],
  voids_per_cluster_hist={int(k):int(c) for k,c in enumerate(np.bincount(g["kc"])) if c},
  dist_over_Rv_median=float(np.median(dRv)),dist_over_Rv_p16_84=[float(x) for x in np.percentile(dRv,[16,84])],
  hamaus_wall_peak_over_Rv=1.09,angular_Rbar_median=float(np.median(Rbar)),angular_Rbar_iso_median=float(np.median(Riso)),
  mass_vs_voidR_slope=float(slope),void_volume_fraction=float(fvol),void_volfrac_target=p["fV"],
  void_core_density=1+p["dv"],cluster_mass_fraction=float(clmass),cluster_mass_fraction_tinker=float(clmass_t),
  net_cell_contrast_rms=float(np.sqrt(np.mean(net**2))),lcdm_sigma_at_cell_median=float(np.median(sig)),
  Reff_median_Mpc=float(np.median(g["Reff"][vi])),Rv_core_median_Mpc=float(np.median(g["Rv"][vi])))
json.dump(st,open(os.path.join(ME,f"stats_check_{Mmin:.0e}.json"),"w"),indent=1); print(json.dumps(st,indent=1))
plt.rcParams.update({"font.size":11})
fig,ax=plt.subplots(3,3,figsize=(15,15)); ax=ax.ravel()
a=ax[0]; b=np.arange(0,kv[vi].max()+2)-0.5; a.hist(kv[vi],b,color="C0",alpha=.7,label="generated")
a.axvline(kv[vi].mean(),color="C0",lw=2,label=f"generated mean {kv[vi].mean():.2f}"); a.axvline(p["kv_mean"],color="k",ls="--",lw=2,label=f"target mean {p['kv_mean']:.1f} (ours)")
a.set_title("(a) clusters bordering\neach void"); a.set_xlabel("clusters per void")
a=ax[1]; b=np.arange(2,g["kc"].max()+2)-0.5; a.hist(g["kc"],b,color="C1",alpha=.7,label="generated")
a.axvspan(4,5,color="green",alpha=.2,label="target 4–5 (from Codis+18 κ≈4 local, ≈6 global)")
a.axvline(g["kc"].mean(),color="C1",lw=2,label=f"generated mean {g['kc'].mean():.2f}"); a.set_title("(b) voids touching\neach cluster"); a.set_xlabel("voids per cluster")
a=ax[2]; a.semilogx(np.sqrt(mb[:-1]*mb[1:]),kcm,"o-",color="C1",label="generated (Aragón-Calvo+10, Codis+18: should rise)")
a.set_title("(b') voids per cluster vs cluster mass\n(literature: rises with mass)"); a.set_xlabel("cluster mass M200m (Msun)"); a.set_ylabel("mean voids per cluster")
a=ax[3]; a.hist(dRv,60,range=(0,3),color="C2",alpha=.7,label="generated cluster–void-centre distance")
a.axvline(1.09,color="k",ls="--",lw=2,label="compensation-wall peak 1.09 Rv (Hamaus+14 fit)"); a.axvline(1,color="gray",lw=1)
a.set_title("(c) cluster distance\nfrom void centre"); a.set_xlabel("distance / Rv (Rv = equal-volume radius, as VIDE)")
a=ax[4]; a.hist(Rbar,50,range=(0,1),color="C3",alpha=.6,density=True,label="generated")
a.hist(Riso,50,range=(0,1),histtype="step",color="k",lw=2,density=True,label="random directions, same counts")
a.set_title("(d) angular spread of a void's clusters\n(0 = evenly around the void, 1 = all on one side)"); a.set_xlabel("resultant length of cluster directions")
a=ax[5]; i=rng.choice(len(g["MC"]),min(8000,len(g["MC"])),replace=False); a.loglog(mR[i],g["MC"][i],".",ms=2,color="C4",alpha=.4,label="generated clusters")
xx=np.logspace(np.log10(mR.min()),np.log10(mR.max()),10); a.loglog(xx,10**np.polyval(np.polyfit(np.log10(mR),np.log10(g["MC"]),1),np.log10(xx)),"k--",label=f"fit slope {slope:.2f} (no literature target; ours)")
a.set_title("(e) cluster mass vs\nsize of its bordering voids"); a.set_xlabel("mean Rv of bordering voids (Mpc)"); a.set_ylabel("M200m (Msun)")
a=ax[6]; bb=np.logspace(np.log10(Mmin),np.log10(g["MC"].max()*1.01),20); h,_=np.histogram(g["MC"],bb); mc=np.sqrt(bb[1:]*bb[:-1])
a.loglog(mc,h/np.diff(np.log(bb))/V,"o",color="C5",label="generated"); mm=np.logspace(np.log10(Mmin),np.log10(g["MC"].max()),50)
a.loglog(mm,lss.cluster_dn_dlnM(mm),"k-",label="Tinker+08 (target)"); a.set_title("(g) cluster mass\nfunction"); a.set_xlabel("M200m (Msun)"); a.set_ylabel("dn/dlnM (Mpc$^{-3}$)")
a=ax[7]; a.hist(g["Reff"][vi],35,color="C6",alpha=.7,label="void Rv (equal-volume radius)"); a.hist(g["Rv"][vi],35,color="C7",alpha=.5,label="non-overlapping core sphere radius")
a.set_title(f"(g) void sizes\nvoid volume fraction {fvol:.2f} (target 0.77, Cautun+14)"); a.set_xlabel("radius (Mpc)")
a=ax[8]; a.hist(net,60,color="C8",alpha=.7,density=True,label="generated net contrast of each void cell")
xs=np.linspace(net.min(),net.max(),200); s0=np.median(sig); a.plot(xs,np.exp(-xs**2/2/s0**2)/np.sqrt(2*np.pi)/s0,"k-",label=f"LCDM σ at cell size ({s0:.3f})")
a.set_title("(f) compensation\nnet contrast of each void + its ridge"); a.set_xlabel("(cell mass − mean)/mean")
for a in ax:
    a.legend(loc="upper center",bbox_to_anchor=(0.5,-0.2),fontsize=9,frameon=False)
fig.suptitle(f"Generated network vs rule targets  (clusters above {Mmin:.0e} Msun, sphere of {R/1000:.0f} Gpc)",fontsize=14)
fig.tight_layout(rect=(0,0,1,0.97),h_pad=2.5); fig.savefig(os.path.join(ME,"stats_check.png" if Mmin==3e15 else f"stats_check_{Mmin:.0e}.png"),dpi=110)
