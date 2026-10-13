import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_vel_collapse_combo_c.py").read().split("sets=")[0])
V=0.1; M=1.0   # held fixed: our infall speed 0.1c, mass 1x7.5e52 kg; D = x*M, R = y*V
X=np.geomspace(10,600,50); Y=np.geomspace(10,500,45)
vesc=lambda Dv: np.sqrt(2*GM0*M/(Dv*GPC))/c
Z=np.full((len(Y),len(X)),np.nan); code=np.zeros_like(Z)
for j,y in enumerate(Y):
    for i,x in enumerate(X):
        Dv=x*M; R=y*V; ve=vesc(Dv)
        if V>=ve: code[j,i]=1; continue
        if R>=0.95*Dv: code[j,i]=2; continue
        rst=Dv/(1-(V/ve)**2); yy=curve(M,Dv,R,rst); Z[j,i]=100*(w*(yy>0)).sum()/w.sum()
print("range %.1f–%.1f"%(np.nanmin(Z),np.nanmax(Z)))
f,ax=plt.subplots(figsize=(7.5,5.8))
ax.pcolormesh(X,Y,np.where(code==1,1,np.nan),cmap=matplotlib.colors.ListedColormap(["0.6"]),shading="auto")
ax.pcolormesh(X,Y,np.where(code==2,1,np.nan),cmap=matplotlib.colors.ListedColormap(["0.85"]),shading="auto")
pc=ax.pcolormesh(X,Y,Z,vmin=0,vmax=100,cmap="RdBu_r",shading="auto")
cs=ax.contour(X,Y,Z,[50],colors="k",linewidths=1.6); ax.clabel(cs,fmt="50%%")
for p in cs.allsegs[0]: print(np.round(p[::max(1,len(p)//8)],1).tolist())
ax.set_xscale("log"); ax.set_yscale("log")
ax.text(0.98,0.04,"v > escape speed\n(masked)",transform=ax.transAxes,ha="right",fontsize=9,bbox=dict(fc="0.6",ec="none"))
ax.text(0.02,0.96,"R ≥ 0.95 D (masked)",transform=ax.transAxes,va="top",fontsize=9,bbox=dict(fc="0.85",ec="none"))
ax.set_xlabel("D / M   (Gpc per 7.5×10⁵² kg)"); ax.set_ylabel("R / |v|   (Gpc per c)")
plt.colorbar(pc,label="% of sky receding (solid-angle weighted)")
ax.set_title("Collapse (free fall from stall radius): % sky receding\nheld fixed: M = 1× (D varies), |v| = 0.1c (R varies)",fontsize=10)
f.savefig("heat_recede_DM_Rv.png",dpi=120,bbox_inches="tight")
