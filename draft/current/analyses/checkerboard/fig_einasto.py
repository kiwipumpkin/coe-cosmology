import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
fig=plt.figure(figsize=(9,15)); gs=fig.add_gridspec(3,1,height_ratios=[1,1,1.25],hspace=0.42)
# --- (a) power spectrum, schematic
k=np.logspace(np.log10(0.01),np.log10(0.4),600)
def bbks(k,G):
    q=k/G; T=np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**-0.25; return k*T**2
cdm=bbks(k,0.5); cdm/=cdm.max()
k0=0.052; bump=1+0.4*np.exp(-0.5*(np.log(k/k0)/0.12)**2)
ax=fig.add_subplot(gs[0]); ax.loglog(k,cdm,"k-",lw=1.6,label="standard CDM shape (h=0.5, Ω=1), as the comparison line in their Fig 2")
ax.loglog(k,cdm*bump,"r-",lw=2.2,label="observed-cluster spectrum, schematic: single peak at k₀ = 0.052 h Mpc⁻¹, ≈1.4× above CDM")
ax.axvline(k0,color="r",ls=":",lw=1); ax.axvspan(2*np.pi/135,2*np.pi/105,color="r",alpha=0.08,label="λ₀ = 2π/k₀ = 120 ± 15 h⁻¹ Mpc")
ax.set_xlabel("wavenumber k  (h Mpc⁻¹)"); ax.set_ylabel("P(k)  (arbitrary units)")
sec=ax.secondary_xaxis("top",functions=(lambda x:2*np.pi/x,lambda x:2*np.pi/x)); sec.set_xlabel("wavelength 2π/k  (h⁻¹ Mpc)")
ax.set_title("(a) Power spectrum: recreation of the shape described for Einasto et al. 1997, Fig 2 (schematic, no data points)",fontsize=10)
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.16),fontsize=8.5,frameon=False)
# --- (b) correlation function, schematic
r=np.linspace(5,400,800); P0=120
xi=(r/20)**-1.8*0.6+0.08*np.exp(-r/350)*np.cos(2*np.pi*r/P0)
ax=fig.add_subplot(gs[1]); ax.plot(r,xi,"b-",lw=2,label="oscillating cluster correlation function (schematic)")
ax.axhline(0,color="k",lw=0.6)
for j in range(1,4): ax.axvline(j*P0,color="b",ls=":",lw=0.8)
ax.annotate("",xy=(240,0.075),xytext=(120,0.075),arrowprops=dict(arrowstyle="<->"))
ax.text(180,0.08,"period = step size = 120 h⁻¹ Mpc",ha="center",fontsize=9)
ax.annotate("",xy=(180,-0.06),xytext=(120,-0.06),arrowprops=dict(arrowstyle="<->"))
ax.text(150,-0.075,"max → min = ½ period = 60",ha="center",va="top",fontsize=9)
ax.set_ylim(-0.12,0.25); ax.set_xlim(0,400); ax.set_xlabel("cluster separation r  (h⁻¹ Mpc)"); ax.set_ylabel("ξ(r)")
ax.set_title("(b) Correlation function: alternating secondary maxima and minima (as described in the text; schematic)",fontsize=10)
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.14),fontsize=8.5,frameon=False)
# --- (c) supercluster-void grid
ax=fig.add_subplot(gs[2]); s=120
for x in range(-360,361,s): ax.axvline(x,color="0.6",lw=0.8); ax.axhline(x,color="0.6",lw=0.8)
for i in range(-3,3):
    for j in range(-3,3): ax.add_patch(Rectangle((i*s+12,j*s+12),s-24,s-24,color="#3b6fd8",alpha=0.12,lw=0))
rng=np.random.default_rng(3)
edges=[((-360,-120),(-240,-120)),((-240,-120),(-240,0)),((-240,0),(-120,0)),((-120,0),(0,0)),((0,0),(0,120)),((0,120),(120,120)),((120,120),(240,120)),((120,-240),(120,-120)),((120,-120),(240,-120)),((0,-240),(120,-240)),((-120,240),(0,240)),((-240,240),(-120,240)),((240,0),(240,120)),((-120,-360),(-120,-240)),((-360,120),(-240,120))]
for (a,b) in edges: ax.plot([a[0],b[0]],[a[1],b[1]],color="#d62728",lw=5,alpha=0.7,solid_capstyle="round")
ax.plot(0,0,"k*",ms=14); ax.text(8,-22,"us",fontsize=10)
ax.set_xlim(-360,360); ax.set_ylim(-360,360); ax.set_aspect("equal")
ax.set_xlabel("supergalactic X  (h⁻¹ Mpc)"); ax.set_ylabel("supergalactic Z  (h⁻¹ Mpc)")
ax.set_title("(c) Grid picture in the style of Einasto et al. 1997 Fig 1: step 120 h⁻¹ Mpc (≈171 Mpc for h=0.7)\nschematic layout only, not the Abell catalogue",fontsize=10)
from matplotlib.lines import Line2D
h=[Line2D([],[],color="#d62728",lw=5,alpha=0.7),Rectangle((0,0),1,1,color="#3b6fd8",alpha=0.15),Line2D([],[],color="0.6")]
ax.legend(h,["chains of rich superclusters along grid lines (illustrative which lines)","voids of almost equal size between them","grid, step 120 ± 20 h⁻¹ Mpc"],loc="upper center",bbox_to_anchor=(0.5,-0.09),fontsize=8.5,frameon=False)
fig.text(0.5,0.01,"Source: Einasto et al. 1997, Nature 385, 139, 'A 120-Mpc periodicity in the three-dimensional distribution of galaxy superclusters', arXiv:astro-ph/9701018.\nAll curves are recreations from the paper's stated numbers (k₀, λ₀, 1.4× excess, slopes n=1 and n=−1.8); no measured points are plotted.",ha="center",fontsize=8)
fig.savefig("fig1_einasto_recreation.png",dpi=140,bbox_inches="tight")
