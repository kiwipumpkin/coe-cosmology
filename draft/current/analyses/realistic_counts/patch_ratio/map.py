import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from ps_core import *
plt.rcParams.update({"font.size":13})
gs=CONFIGS["in10"](); g=gs[0]; z=np.load("thr_in10.npz"); d=z[z.files[0]]; P=g.P; T=json.load(open("thr_in10.json")); B=T.get("B_1x_primal",T["B_star"])
p=pinfo("extreme 50 Mpc"); f0=p["f0"]; R0=f0/(1-f0); f=f0+d/p["Delta"]; rel=(f/(1-f))/R0
lon=np.degrees(np.arctan2(P[:,1],P[:,0]))%360; r=np.linalg.norm(P,axis=1)
# --- coherence: azimuthal m-content of the ring-plane field (r>9 Gpc, |z|<1) and 3D correlation length
sel=(np.abs(P[:,2])<1)&(r>9); bins=np.arange(0,361,5); prof=np.array([d[sel&(lon>=a)&(lon<a+5)].mean() for a in bins[:-1]])
pw=np.abs(np.fft.rfft(prof-prof.mean()))**2; pw=pw/pw[1:].sum(); mfrac={int(m):float(pw[m]) for m in range(1,11)}
F=np.zeros(g.mask.shape); F[g.mask]=d; W=g.mask.astype(float); n=F.shape[0]
c=np.fft.irfftn(np.abs(np.fft.rfftn(F,s=(2*n,)*3))**2,s=(2*n,)*3); cw=np.fft.irfftn(np.abs(np.fft.rfftn(W,s=(2*n,)*3))**2,s=(2*n,)*3)
ii=np.fft.fftfreq(2*n,1/(2*n))*g.H; RR=np.sqrt(ii[:,None,None]**2+ii[None,:,None]**2+ii[None,None,:]**2)
xi=np.where(cw>50,c/np.maximum(cw,1e-9),np.nan); rb=np.arange(0,14,0.5); xr=np.array([np.nanmean(xi[(RR>=a)&(RR<a+0.5)]) for a in rb]); xr/=xr[0]
half=float(rb[np.argmax(xr<0.5)]); zero=float(rb[np.argmax(xr<0)])
# --- lobes
lobes=[]
for cc in [15,75,128,162,195,228,263,315]:
    dl=np.abs((lon-cc+180)%360-180); k=(dl<22)&(r>9)
    sg=np.sign(d[k].mean()); kk=k&(sg*d>0.1*B)
    if not kk.any(): continue
    lobes.append(dict(lon=cc,sign=int(sg),rel_ratio=float(rel[kk].mean()),df=float((f[kk]-f0).mean()),mean_delta=float(d[kk].mean()),
                      vol=float(kk.sum()*g.H**3),frac_cells_at_bound=float((np.abs(d[kk])>0.98*B).mean())))
json.dump(dict(B=B,lobes=lobes,m_power_frac=mfrac,corr_half_gpc=half,corr_zero_gpc=zero,xi=dict(zip(rb.tolist(),xr.tolist())),
               vol_active=float((np.abs(d)>0.1*B).sum()*g.H**3),vol_at_bound=float((np.abs(d)>0.98*B).sum()*g.H**3)),open("pattern_in10.json","w"),indent=1)
# --- map
s=(P[:,2]>0)&(P[:,2]<0.6); fig,ax=plt.subplots(figsize=(10.5,10))
v=100*(rel-1); lim=np.ceil(np.abs(v[s]).max())
sc=ax.scatter(P[s,0],P[s,1],c=v[s],cmap="RdBu_r",vmin=-lim,vmax=lim,s=240,marker="s",linewidths=0)
ax.add_patch(Circle((0,0),14.3,fill=False,ls="--",color="k",lw=1.5)); ax.plot(0,0,"k+",ms=16,mew=2.5)
for l in lobes:
    a=np.radians(l["lon"]); ax.text(16.3*np.cos(a),16.3*np.sin(a),f"over:under\n{l['rel_ratio']:.2f} : 1",ha="center",va="center",fontsize=12,
        color="#922b21" if l["sign"]>0 else "#1a5276",bbox=dict(facecolor="white",alpha=0.85,edgecolor="none",pad=1.5))
ax.set_xlim(-19.5,19.5); ax.set_ylim(-19.5,19.5); ax.set_aspect("equal"); ax.set_xlabel("Gpc"); ax.set_ylabel("Gpc")
cb=fig.colorbar(sc,ax=ax,shrink=0.7,pad=0.02); cb.set_label("over : under patch ratio vs cosmic mean (%)")
ax.set_title(f"Required over:under patch ratio (relative to cosmic mean), ring-plane slice\ninside, 1 Gpc domains, 50 Mpc patches (+0.83 / −0.39)\ndashed = ring; labels = lobe means in the outer shell (r > 9 Gpc)",fontsize=13)
fig.tight_layout(); fig.savefig("map_over_under_ratio.png",dpi=110); plt.close(fig)
print(json.dumps(dict(lobes=lobes,m=mfrac,half=half,zero=zero),indent=0))
