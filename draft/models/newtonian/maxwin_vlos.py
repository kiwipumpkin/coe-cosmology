import numpy as np, matplotlib.pyplot as plt
exec(open('maxwin.py').read().split('# --- family A')[0])
plt.rcParams.update({"font.size":17,"axes.labelsize":18,"axes.titlesize":18,"xtick.labelsize":16,"ytick.labelsize":16,"legend.fontsize":14})
def go(rs,RD,nm):
  r,v,_=fieldB(rs); F=lambda x:np.where(x<r[-1],np.interp(x,r,v),np.nan); Rs=np.linspace(0,RD,300); xo=np.array([1.,0])
  f,ax=plt.subplots(figsize=(8,6))
  for a,cl in zip([0,45,90,135,180],["#b2182b","#ef8a62","#4d4d4d","#67a9cf","#2166ac"]):
    nn=np.array([np.cos(np.radians(a)),np.sin(np.radians(a))]); xp=xo+np.outer(Rs,nn); rp=np.linalg.norm(xp,axis=1)
    vl=((F(rp)[:,None]*xp/rp[:,None]-F(1)*xo)*nn).sum(1)*vu; ax.plot(Rs*3200,vl,color=cl,lw=3.5,label=f"{a}°")
    print(nm,a,np.nanmax(vl),np.nanmin(vl))
  ax.plot(Rs*3200,-Rs*vu,"k--",lw=2.5,label="no mass"); ax.axhline(0,color="gray",lw=1)
  ax.set_xlim(0,RD*3200); ax.set_xlabel("galaxy distance R (Mpc)"); ax.set_ylabel("v_LOS (units of c)")
  ax.set_title(f"Collapsing, stall/D = {rs:.2f}"); ax.legend(frameon=False,ncol=2,fontsize=13); f.tight_layout(); f.savefig(f"maxwin_vlos_{nm}.png",dpi=110)
go(1.05,0.049,"stall105"); go(1.35,0.246,"stall135"); go(2.5,0.3,"stall250_vc")
B=np.load('maxwin_B.npy',allow_pickle=True).item(); j=np.argmin(abs(B['RDb']-.3)); i=np.argmin(abs(B['rss']-2.5)); print("win",B['WB'][i,j])
o=evalf(*fieldB(2.5),0.3); print(o['vmax'],o['max']['win'])
