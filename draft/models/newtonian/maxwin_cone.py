import numpy as np, matplotlib.pyplot as plt
exec(open('maxwin.py').read().split('# --- family A')[0])
plt.rcParams.update({"font.size":17,"axes.labelsize":18,"axes.titlesize":18,"xtick.labelsize":16,"ytick.labelsize":16,"legend.fontsize":13})
th=np.radians(np.linspace(0.05,179.95,720)); ph=np.radians(np.linspace(0,360,720,endpoint=False)); TH,PH=np.meshgrid(th,ph)
d=np.stack([np.cos(TH),np.sin(TH)*np.cos(PH),np.sin(TH)*np.sin(PH)],-1); W=np.sin(TH)
cens=[0,45,90,135,180]; masks=[np.degrees(np.arccos(np.clip(d@np.array([np.cos(np.radians(a)),np.sin(np.radians(a)),0]),-1,1)))<75 for a in cens]
nn=np.stack([np.cos(th),np.sin(th)],1); cols=["#b2182b","#ef8a62","#4d4d4d","#67a9cf","#2166ac"]
for rs,RD,nm in [(1.05,.049,"stall105"),(1.35,.246,"stall135"),(2.5,.3,"stall250_vc")]:
  r,v,_=fieldB(rs); F=lambda x:np.where(x<r[-1],np.interp(x,r,v),np.nan); Rs=np.linspace(0,RD,120); X=Rs*3200
  out=np.zeros((6,len(Rs))); Rinv=None
  for k,R in enumerate(Rs):
    xp=np.array([1.,0])+R*nn; rp=np.linalg.norm(xp,axis=1); vl=((F(rp)[:,None]*xp/rp[:,None]-F(1)*np.array([1.,0]))*nn).sum(1)*vu
    V=np.broadcast_to(vl,TH.shape); ok=np.isfinite(V)
    if Rinv is None and not ok.all(): Rinv=X[k]
    for i,mk in enumerate(masks+[np.ones_like(W,bool)]):
      g=mk&ok; out[i,k]=(W*np.where(g,V,0)).sum()/(W*g).sum() if g.any() else np.nan
  lim=X<(Rinv if Rinv else 1e9); Hs=[np.polyfit(X[lim&np.isfinite(o)],o[lim&np.isfinite(o)]*c,1)[0] for o in out]
  print(nm,"fit R<%s Mpc |"%(f"{Rinv:.0f}" if Rinv else "all")," ".join(f"{a}°:{h:.0f}" for a,h in zip(cens+["sky"],Hs)))
  f,ax=plt.subplots(figsize=(8,6))
  for i,a in enumerate(cens): ax.plot(X,out[i],lw=3.5,color=cols[i],label=f"cone {a}° (H≈{Hs[i]:.0f})")
  ax.plot(X,out[5],lw=4,color="#1b7837",label=f"whole sky (H≈{Hs[5]:.0f})"); ax.plot(X,-X*70/c,"k--",lw=2.5,label="no mass (−70)")
  if Rinv: ax.axvline(Rinv,color="m",ls=":",lw=2.5)
  ax.axhline(0,color="gray",lw=1); ax.set_xlim(0,X[-1]); ax.set_xlabel("galaxy distance R (Mpc)"); ax.set_ylabel("cone-average v_LOS (c)")
  ax.set_title(f"Collapsing, stall/D = {rs:.2f}, 75° cones"); ax.legend(frameon=False,fontsize=12,ncol=2,loc="upper left"); f.tight_layout(); f.savefig(f"maxwin_cone75_{nm}.png",dpi=110); plt.close(f)
