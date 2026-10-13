import numpy as np, matplotlib.pyplot as plt
exec(open('maxwin.py').read().split('# --- family A')[0])
plt.rcParams.update({"font.size":17,"axes.labelsize":18,"axes.titlesize":18,"xtick.labelsize":16,"ytick.labelsize":16,"legend.fontsize":14})
B=np.load('maxwin_B.npy',allow_pickle=True).item(); rss,RDb,WB=B['rss'],B['RDb'],B['WB']
VB=np.full_like(WB,np.nan); WBall=np.full_like(WB,np.nan)
for i,rs in enumerate(rss):
  fb=fieldB(rs)
  for j,RD in enumerate(RDb):
    o=evalf(*fb,min(RD,rs-1.0101)) if rs>1+RD+.02 else None
    if o: VB[i,j]=o['vmax']
def best(W,V,lim):
  Wc=np.where(V<lim,W,np.nan); return np.unravel_index(np.nanargmax(Wc),W.shape)
def rep(tag,fa,RD,extra):
  o=evalf(*fa,RD); m=o['max']; s=o['small']
  print(f"{tag}: {extra} R/D={RD:.3f} | Rmax win {m['win']:.1f}% (plain {m['plain']:.1f}%) mean {m['mean']*70:.1f} sd {m['sd']*70:.1f} H0/90/180 {m['y'][0]*70:.0f}/{m['y'][360]*70:.0f}/{m['y'][720]*70:.0f} | smallR win {s['win']:.1f}% plain {s['plain']:.1f}% | vmax {o['vmax']:.2f}c (D=3.2Gpc)")
  return o
cases=[]
i,j=best(WB,VB,1e9); cases.append(("stall, no c-limit",fieldB(rss[i]),RDb[j],f"rs/D={rss[i]:.2f}"))
i,j=best(WB,VB,1.0); cases.append(("stall, v<c",fieldB(rss[i]),RDb[j],f"rs/D={rss[i]:.2f}"))
cases.append(("stall, R/D=0.246",fieldB(rss[np.nanargmax(WB[:,np.argmin(abs(RDb-.246))])]),RDb[np.argmin(abs(RDb-.246))],f"rs/D={rss[np.nanargmax(WB[:,np.argmin(abs(RDb-.246))])]:.2f}"))
A=np.load('maxwin_A2.npy',allow_pickle=True).item(); W,V=A['res'][0.246]; i,j=best(W,V,1e9)
cases.append(("start v=-kr, no c-limit",fieldA(A['mus'][j],A['taus'][i]),0.246,f"M={A['mus'][j]/mu1:.0f}x kT={A['taus'][i]:.2f}"))
outs=[rep(*c) for c in cases]
for nm,(c,o) in zip(["stall_free","stall_vc","stall_RD0246","startv_free"],zip(cases,outs)):
  f,ax=plt.subplots(figsize=(8,6)); ax.plot(ang,o['max']['y']*70,lw=3.5,color="#2166ac",label="best"); ax.axhline(-70,color="k",ls="--",lw=2.5,label="no mass")
  ax.axhline(0,color="gray",lw=1); ax.set_xlim(0,360); ax.set_xticks(range(0,361,90)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_ylabel("dv/dR (km/s/Mpc)")
  ax.set_title(f"Collapsing: {c[0]}\n{c[3]}, R/D={c[2]:.2f}, window {o['max']['win']:.0f}%",fontsize=16); ax.legend(frameon=False); f.tight_layout(); f.savefig(f"maxwin_slope_{nm}.png",dpi=110); plt.close(f)
# contours
f,ax=plt.subplots(figsize=(8,6)); Wf=np.nan_to_num(WB,nan=0)
cs=ax.contourf(RDb,rss,WB,levels=np.arange(0,75,5),cmap="viridis"); f.colorbar(cs,label="recession window (% sky)")
ax.fill_between(RDb,1,1+RDb+.02,color="0.6"); ax.contour(RDb,rss,VB,levels=[1],colors="w",linestyles="--",linewidths=3)
ax.text(.5,1.2,"shells outside stall\n(not all inward)",fontsize=13,ha="center")
ax.set_xlabel("R/D (view / distance to CoE)"); ax.set_ylabel("stall radius / D"); ax.set_title("Stall collapse: window at edge of view\nwhite dashed: v = c (D = 3.2 Gpc)",fontsize=16)
f.tight_layout(); f.savefig("maxwin_contour_stall.png",dpi=110); plt.close(f)
f,ax=plt.subplots(figsize=(8,6)); cs=ax.contourf(A['mus']/mu1,A['taus'],W,levels=np.arange(0,25,2),cmap="viridis"); f.colorbar(cs,label="recession window (% sky)")
ax.contour(A['mus']/mu1,A['taus'],V,levels=[1],colors="w",linestyles="--",linewidths=3)
ax.set_xlabel("mass (× 1×)"); ax.set_ylabel("kT (collapse time × k)"); ax.set_title("Start v=−kr collapse, R/D=0.246\nwhite dashed: v = c (D = 3.2 Gpc)",fontsize=16)
f.tight_layout(); f.savefig("maxwin_contour_startv.png",dpi=110)
