import numpy as np, matplotlib.pyplot as plt
exec(open('maxwin.py').read().split('# --- family A')[0])
plt.rcParams.update({"font.size":17,"axes.labelsize":18,"axes.titlesize":18,"xtick.labelsize":16,"ytick.labelsize":16,"legend.fontsize":14})
xo=np.array([1.,0]); res=[]
for rs,RD,nm in [(1.05,.049,"stall105"),(1.35,.246,"stall135"),(2.5,.3,"stall250_vc")]:
  r,v,_=fieldB(rs); F=lambda x:np.where(x<r[-1],np.interp(x,r,v),np.nan); Rs=np.linspace(0,RD,200); m=[];bad=[]
  for R in Rs:
    xp=xo+R*n; rp=np.linalg.norm(xp,axis=1); vl=((F(rp)[:,None]*xp/rp[:,None]-F(1)*xo)*n).sum(1)*vu
    ok=np.isfinite(vl); m.append((w*np.where(ok,vl,0)).sum()/w[ok].sum()); bad.append(100*(w*~ok).sum()/w.sum())
  m=np.array(m); bad=np.array(bad); Mpc_=Rs*3200; H=np.polyfit(Mpc_,m*c,1)[0]; Hv=np.polyfit(Mpc_[bad==0],m[bad==0]*c,1)[0] if (bad==0).sum()>2 else np.nan
  Rfirst=Mpc_[np.argmax(bad>0)] if bad.max()>0 else None
  print(nm,"H fit %.1f, H fit valid-only-R %.1f, end v %.3fc, first invalid R %s, max invalid %.1f%%"%(H,Hv,m[-1],Rfirst,bad.max()))
  res.append((rs,Mpc_,m,H,Rfirst))
  f,ax=plt.subplots(figsize=(8,6)); ax.plot(Mpc_,m,lw=3.5,color="#2166ac",label="sky average (H≈%.0f)"%H); ax.plot(Mpc_,-Mpc_*70/c,"k--",lw=2.5,label="no mass (−70)")
  if Rfirst: ax.axvline(Rfirst,color="#b2182b",ls=":",lw=2.5,label="beyond: some sky past stall")
  ax.axhline(0,color="gray",lw=1); ax.set_xlim(0,Mpc_[-1]); ax.set_xlabel("galaxy distance R (Mpc)"); ax.set_ylabel("sky-average v_LOS (c)")
  ax.set_title(f"Collapsing, stall/D = {rs:.2f}"); ax.legend(frameon=False,fontsize=13); f.tight_layout(); f.savefig(f"maxwin_skyavg_{nm}.png",dpi=110); plt.close(f)
f,ax=plt.subplots(figsize=(8,6))
for (rs,X,m,H,Rf),cl in zip(res,["#b2182b","#2166ac","#1b7837"]): ax.plot(X,m,lw=3.5,color=cl,label=f"stall/D={rs:.2f} (H≈{H:.0f})")
X=np.linspace(0,960,10); ax.plot(X,-X*70/c,"k--",lw=2.5,label="no mass (−70)"); ax.axhline(0,color="gray",lw=1)
ax.set_xlim(0,960); ax.set_xlabel("galaxy distance R (Mpc)"); ax.set_ylabel("sky-average v_LOS (c)"); ax.set_title("Collapsing, sky average"); ax.legend(frameon=False,fontsize=13)
f.tight_layout(); f.savefig("maxwin_skyavg_all.png",dpi=110)
