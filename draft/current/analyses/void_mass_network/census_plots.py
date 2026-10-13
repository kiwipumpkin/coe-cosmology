"""Plain census charts from results/census_arrays.npz + census_3e+15.json: after the non-crossing rule = filled,
before it = black step outline. Legends below."""
import os,json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ME=os.path.dirname(os.path.abspath(__file__)); a=np.load(os.path.join(ME,"results","census_arrays.npz")); c=json.load(open(os.path.join(ME,"census_3e+15.json")))
plt.rcParams.update({"font.size":11}); H=0.7
def leg(ax): ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.2),fontsize=9.5,frameon=False)
def ih(x,key,col,lab,lo=0):
    A=a["after_"+key]; B=a["before_"+key]; b=np.arange(lo,max(A.max(),B.max())+2)-0.5
    x.hist(A,b,color=col,alpha=.65,label=f"ours, non-crossing: mean {A.mean():.2f}, sd {A.std():.2f}")
    x.hist(B,b,histtype="step",color="k",lw=1.3,label=f"ours, before the rule: mean {B.mean():.2f}, sd {B.std():.2f}"); x.set_ylabel("number of nodes")
fig,ax=plt.subplots(3,2,figsize=(13,16.5)); ax=ax.ravel()
x=ax[0]; ih(x,"deg_ccall","#c0392b","filaments per cluster (direct + via one empty node)")
x.axvspan(1,4,color="0.5",alpha=.15,label="Colberg+05 (astro-ph/0406665): most clusters have 1–4 filaments, <10% none")
x.axvline(2,color="k",ls="--",lw=1.5,label="Aragón-Calvo+10 (1007.0742): ≈2 filaments at 10¹⁴ Msun/h")
x.axvline(5,color="k",ls=":",lw=2,label="Aragón-Calvo+10: ≈5 at ≥10¹⁵ Msun/h; Colberg+05: almost 5 for the most massive")
x.set_title("cluster–cluster (filaments), all\nper cluster"); x.set_xlabel("filaments per cluster"); leg(x)
x=ax[1]; ih(x,"deg_cc","#e74c3c","direct foam-edge filaments only")
x.axvline(4,color="k",ls="--",lw=1.5,label="Codis+18 (1803.11477): ≈4 per peak locally in 3D\n(≈6.1 globally; counted over all peaks, not only rare ones)")
x.set_title("cluster–cluster, direct foam edges only\nper cluster (our clusters are 42% of web nodes)"); x.set_xlabel("filaments per cluster"); leg(x)
x=ax[2]; ih(x,"deg_cv_c","#8e44ad","",lo=0); x.axvspan(4,5,color="green",alpha=.18,label="derived from Codis+18: 2+κ/2 = 4–5")
x.set_title("cluster–void, per cluster"); x.set_xlabel("voids per cluster"); leg(x)
x=ax[3]; ih(x,"deg_cv_v","#2980b9",""); x.set_title("cluster–void, per void\n(no published count found)"); x.set_xlabel("clusters per void"); leg(x)
x=ax[4]; ih(x,"deg_vv","#16a085","",lo=0); x.axvline(15.54,color="k",ls="--",lw=1.5,label="Poisson–Voronoi foam (exact): 48π²/35+2 = 15.54 neighbours")
x.set_title("void–void (shared walls), per void"); x.set_xlabel("neighbour voids per void"); leg(x)
x=ax[5]
for nm,key,bk,col,mk in (("filaments per cluster vs cluster-mass third","deg_ccall","mbin","#c0392b","o-"),("voids per cluster vs cluster-mass third","deg_cv_c","mbin","#8e44ad","o-"),
                        ("clusters per void vs void-size third","deg_cv_v","vbin","#2980b9","s--"),("neighbour voids per void vs void-size third","deg_vv","vbin","#16a085","s--")):
    A=a["after_"+key]; b=a[bk]; x.plot([0,1,2],[A[b==k].mean() for k in range(3)],mk,color=col,label="ours, non-crossing: "+nm)
x.set_xticks([0,1,2]); x.set_xticklabels(["smallest third","middle","largest third"]); x.set_ylabel("mean connections"); x.set_ylim(0,None)
x.set_title("split by cluster mass / void size\n(literature: filaments per cluster rise with mass)"); leg(x)
fig.suptitle("Connections per node, all three link types (one network, clusters >3e15 Msun; nodes within 4.5 Gpc)",fontsize=14)
fig.tight_layout(rect=(0,0,1,.975),h_pad=2); fig.savefig(os.path.join(ME,"census_counts.png"),dpi=100); plt.close(fig)
# lengths
T=[("cc","#e74c3c","cluster–cluster, direct"),("cc2","#e67e22","cluster–cluster, via one empty node"),("cv","#8e44ad","cluster–void"),("vv","#16a085","void–void")]
fig,ax=plt.subplots(2,1,figsize=(11,15))
x=ax[0]; bins=np.logspace(0.5,3.5,70)
for k,col,lab in T:
    A=a[f"after_L_{k}"]; B=a[f"before_L_{k}"]; x.hist(A,bins,color=col,alpha=.35,label=f"ours, non-crossing, {lab}: median {np.median(A):.0f} Mpc")
    x.hist(B,bins,histtype="step",color=col,lw=1.6,label=f"   before the rule: median {np.median(B):.0f} Mpc")
x.axvspan(5,6,color="0.6",alpha=.35,label="Galárraga-Espinosa+20 (2003.09697): filament lengths peak at ≈5–6 Mpc (TNG300-1)")
x.axvspan(5/H,15/H,color="#c0392b",alpha=.12,label="Colberg+05: 2/3 of straight cluster–cluster filaments are 5–15 Mpc/h")
x.axvline(65.6,color="k",ls="--",lw=1.5,label="Galárraga-Espinosa+20: longest TNG300-1 filament 65.6 Mpc")
x.set_xscale("log"); x.set_xlabel("link length (Mpc)"); x.set_ylabel("number of links")
x.set_title("Link lengths in Mpc\n(our clusters are the rarest, >3e15 Msun, so links are far longer than catalogued filaments)"); leg(x)
x=ax[1]; bins=np.linspace(0,4,80)
for k,col,lab in T:
    A=a[f"after_LR_{k}"]; A=A[np.isfinite(A)]; B=a[f"before_LR_{k}"]; B=B[np.isfinite(B)]
    x.hist(A,bins,color=col,alpha=.35,label=f"ours, non-crossing, {lab}: median {np.median(A):.2f} Rv"); x.hist(B,bins,histtype="step",color=col,lw=1.6,label=f"   before the rule: median {np.median(B):.2f} Rv")
x.axvline(1.09,color="k",ls="--",lw=1.5,label="void ridge peak 1.09 Rv (from the Hamaus+14 profile, 1403.5499)")
x.set_xlabel("link length / Rv  (Rv = equal-volume void radius of the void(s) the link touches or shares)"); x.set_ylabel("number of links"); x.set_title("Link lengths in units of void radius"); leg(x)
fig.tight_layout(h_pad=3); fig.savefig(os.path.join(ME,"census_lengths.png"),dpi=100); plt.close(fig)
# angles
fig,ax=plt.subplots(2,2,figsize=(13,12)); ax=ax.ravel(); th=np.linspace(0,180,200)
for x,k,col,t in ((ax[0],"ang_cc","#c0392b","between filaments at a cluster"),(ax[1],"ang_cvc","#8e44ad","between a cluster's void links"),
                  (ax[2],"ang_vv","#16a085","between a void's wall links"),(ax[3],"ang_cvv","#2980b9","between a void's cluster links")):
    A=a["after_"+k]; B=a["before_"+k]
    x.hist(A,60,range=(0,180),density=True,color=col,alpha=.6,label=f"ours, non-crossing: median {np.median(A):.0f}°")
    x.hist(B,60,range=(0,180),density=True,histtype="step",color="k",lw=1.2,label=f"ours, before the rule: median {np.median(B):.0f}°")
    x.plot(th,np.sin(np.radians(th))/2*np.pi/180,color="0.4",ls="--",lw=1.5,label="random directions (reference)")
    x.axvline(109.47,color="0.3",ls=":",lw=1.5,label="tetrahedral 109.5° (4 links evenly spread)")
    x.set_title(t+"\n(no published distribution found)"); x.set_xlabel("angle between two links at the same node (deg)"); leg(x)
fig.tight_layout(h_pad=2.5); fig.savefig(os.path.join(ME,"census_angles.png"),dpi=100)
