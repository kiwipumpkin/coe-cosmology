import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,"/workspace/clump_scan/inside_visible/smooth_field")
from lcdm import sigma_TH
R=np.logspace(0,np.log10(20000),200)            # Mpc
sig=np.array([sigma_TH(r/1000) for r in R])
np.savetxt("lcdm_sigma_tophat.txt",np.c_[R,sig],header="R_Mpc sigma_tophat_linear_z0 (EH98 no-wiggle, s8=0.81, Om=0.31, h=0.70)")

fig,(a,b)=plt.subplots(2,1,figsize=(11,10),sharex=True,gridspec_kw=dict(hspace=0.06))
for ax in (a,b):
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1,20000)
    ax.fill_between(R,1e-6,3*sig,color="#b9d7ee",zorder=0,label="Standard cosmology: typical range (3σ)")
    ax.fill_between(R,1e-6,sig,color="#6aa6d6",zorder=0,label="Standard cosmology: 1σ")
    ax.axvspan(1,10,color="0.9",zorder=0.5)
    ax.axvline(95,color="0.5",ls=":",lw=1.2)
    ax.grid(True,which="major",color="0.85",lw=0.6)
a.set_ylim(1e-4,3e3); b.set_ylim(3e-5 if False else 1e-4,1.0); b.invert_yaxis()
a.set_ylabel("over-density  δ  (+)",fontsize=13)
b.set_ylabel("under-density  δ  (−)",fontsize=13)
b.set_yticks([1e-4,1e-3,1e-2,0.1,1]); b.set_yticklabels(["−0.0001","−0.001","−0.01","−0.1","−1 (empty)"])
a.set_yticks([1e-4,1e-3,1e-2,0.1,1,10,100,1000]); a.set_yticklabels(["0.0001","0.001","0.01","0.1","1","10","100","1000"])
b.set_xlabel("region radius (Mpc)        1 Gpc = 1000 Mpc",fontsize=13)
b.set_xticks([1,10,100,1000,10000]); b.set_xticklabels(["1","10","100","1,000","10,000"])
kw=dict(fontsize=10.5)
def obs(ax,x,y,lab,dx=1.15,dy=1.0,ha="left",xerr=None,yerr=None):
    ax.errorbar(x,y,xerr=xerr,yerr=yerr,fmt="o",color="k",ms=8,capsize=3,zorder=5)
    ax.annotate(lab,(x,y),xytext=(x*dx,y*dy),ha=ha,va="center",**kw,zorder=6)
# --- observed over-dense
obs(a,2.84,645,"Coma cluster\n(R200)",dx=1.25)
obs(a,11.4,20,"Shapley core",dx=1.25)
obs(a,36.6,1.62,"Saraswati",dx=0.85,dy=2.0,ha="right")
obs(a,79,1.0,"Sloan Great Wall",dx=1.25,dy=1.0)
# --- observed under-dense
obs(b,44,0.75,"Boötes void",dx=0.8,ha="right",yerr=[[0.08],[0.08]])
obs(b,40,0.4,"SDSS voids (lensing avg)",dx=0.5,dy=0.6,ha="right",xerr=[[19],[39]])
obs(b,314,0.14,"Eridanus supervoid (CMB Cold Spot)",dx=1.2,dy=1.25,yerr=[[0.04],[0.04]])
b.errorbar(250,0.24,yerr=[[0.09],[0.09]],fmt="o",mfc="white",color="k",ms=8,capsize=3,zorder=5)
b.annotate("Local Hole (disputed)",(250,0.24),xytext=(330,0.30),ha="left",va="center",**kw)
# --- quasar dipole context
a.plot(7000,0.0154,"D",mfc="white",color="k",ms=8,zorder=5)
a.annotate("Quasar count dipole 1.5%\n(sky tilt in counts,\nnot a density)",(7000,0.0154),xytext=(7000,0.0032),ha="center",va="center",fontsize=9.5,color="0.25")
# --- size-only giant structures
a.annotate("",xy=(200,7),xytext=(1500,7),arrowprops=dict(arrowstyle="<->",color="0.45"))
a.text(550,2.6,"Giant Arc, Big Ring, Huge-LQG, Her–CrB Wall:\nsizes only, no δ measured (debated)",ha="center",va="center",fontsize=9.5,color="0.35")
# --- Chris's required field
gx=np.array([1000,2000,4000]); pk=np.array([0.068,0.055,0.059]); rm=np.array([0.0105,0.012,0.014])
for ax in (a,b):
    ax.plot(gx,pk,"-*",color="#d62728",ms=16,lw=2,zorder=7,label="Required smooth field: peak")
    ax.plot(gx,rm,"--o",color="#d62728",mfc="white",ms=8,lw=1.5,zorder=7,label="Required smooth field: rms")
a.annotate("Required field\npeak ±5–7%",(4000,0.059),xytext=(800,0.16),ha="center",fontsize=11,color="#d62728",weight="bold")
a.annotate("rms ≈1%",(1000,0.0105),xytext=(850,0.0105),ha="right",va="center",fontsize=10.5,color="#d62728")
# --- point-mass fits
a.plot(830,100,"s",color="#7b3fa0",ms=11,zorder=7); a.annotate("Point-mass fit cores\n(+100 over 0.8 Gpc)",(830,100),xytext=(1150,100),va="center",fontsize=10.5,color="#7b3fa0",weight="bold")
b.errorbar(3750,0.7,xerr=[[250],[250]],fmt="s",color="#7b3fa0",ms=11,capsize=3,zorder=7)
b.annotate("Point-mass fit voids\n(−0.7 over 3.5–4 Gpc)",(3750,0.7),xytext=(2900,0.72),ha="right",va="center",fontsize=10.5,color="#7b3fa0",weight="bold")
# --- region notes
a.text(1.25,0.0002,"linear theory\nunder-predicts\nhere (clusters\nare nonlinear)",fontsize=9,color="0.4",va="bottom")
a.text(100,0.00013,"← lumpy | smooth on average →\n   (homogeneity scale ~95 Mpc)",fontsize=9,color="0.4",ha="left")
b.text(1.4,0.00025,"blue = what standard cosmology expects\nfor a region that size today",fontsize=9.5,color="#2c6ea3")
h,l=a.get_legend_handles_labels(); a.legend(h[:2],l[:2],loc="upper right",bbox_to_anchor=(1.0,0.99),fontsize=9.5,framealpha=0.9)
a.set_title("Biggest observed over- and under-densities vs region size",fontsize=15)
fig.savefig("density_vs_scale.png",dpi=130,bbox_inches="tight")
for r in [10,100,300,1000,2000,4000]: print(r, sigma_TH(r/1000), 3*sigma_TH(r/1000))
