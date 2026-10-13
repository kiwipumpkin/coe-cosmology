import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src=open("rel_vel_collapse_combo_c.py").read()
exec(src.split("f,ax=")[0])          # defines curve, sets, w, ang (0..360 step .125)
# axisymmetric about CoE axis: v depends on polar angle th in [0,180]; use 0..180 half of curve
th=ang[ang<=180]; i180=len(th)
tb=np.radians(th); dth=np.gradient(tb)
def kernel(alpha):
    a=np.radians(alpha); K=np.zeros((len(th),len(th)))
    for i,t0 in enumerate(tb):
        s=np.sin(t0)*np.sin(tb); num=np.cos(a)-np.cos(t0)*np.cos(tb)
        with np.errstate(divide="ignore",invalid="ignore"):
            x=np.where(s>1e-12,num/s,np.where(num<=0,-np.inf,np.inf))
        frac=np.arccos(np.clip(x,-1,1))/np.pi       # fraction of ring within cone
        K[i]=frac*np.sin(tb)*dth
    return K/K.sum(1,keepdims=True)
K60,K75=kernel(60),kernel(75)
def smooth(v,K):
    h=K@v[:i180]; return np.concatenate([h,h[::-1][1:len(v)-i180+1]])
wh=np.sin(tb)*dth; wh/=wh.sum()
def stats(vh): return 100*(wh*(vh>0)).sum(), (wh*vh).sum()
names=[]
f,axs=plt.subplots(2,4,figsize=(16,7.5),sharex=True); axs=axs.ravel()
for k,p in enumerate(sets):
    y=curve(**p); y60=smooth(y,K60); y75=smooth(y,K75)
    fr=[stats(z[:i180])[0] for z in (y,y60,y75)]; mn=stats(y[:i180])[1]
    rs="∞" if np.isinf(p["rst"]) else "%.4g"%p["rst"]
    ax=axs[k]; ax.axhline(0,color="k",lw=.8)
    ax.plot(ang,y,"k",lw=2,label="raw (rec %.0f%%)"%fr[0]); ax.plot(ang,y60,"#c0392b",lw=1.8,label="60° cone (rec %.0f%%)"%fr[1]); ax.plot(ang,y75,"#2e5fa8",lw=1.6,ls="--",label="75° cone (rec %.0f%%)"%fr[2])
    ax.set_title(("BASELINE: " if k==0 else "")+"M=%g×, D=%g, R=%g, stall=%s Gpc\nsky mean %+.3g c"%(p["m"],p["Dv"],p["sv"],rs,mn),fontsize=9)
    ax.set_xticks(range(0,361,45)); ax.set_xlim(0,360); ax.tick_params(labelsize=7); ax.legend(fontsize=7,frameon=False)
    print("%s raw %.1f  c60 %.1f  c75 %.1f  mean %.4f"%(p,fr[0],fr[1],fr[2],mn))
for ax in axs[4:]: ax.set_xlabel("angle θ (deg): 0° = away from CoE, 180° = toward",fontsize=8)
for ax in axs[::4]: ax.set_ylabel("v_LOS / c (+ = receding)",fontsize=8)
f.suptitle("Collapse (stall-radius free fall): raw vs Migkas-style cone-averaged v_LOS (solid-angle weighted top-hat cones, 60° and 75°)",fontsize=11)
f.tight_layout(); f.savefig("cone_curves.png",dpi=110)
# heat maps
D=58.2; X=np.linspace(0.02,0.9,45); Y=np.linspace(0,1,41)
Z={k:np.zeros((len(Y),len(X))) for k in ("raw","60","75")}
for j,yv in enumerate(Y):
    for i,x in enumerate(X):
        rst=np.inf if yv>=1 else D/(1-yv**2); v=curve(1,D,x*D,rst)[:i180]
        Z["raw"][j,i]=stats(v)[0]; Z["60"][j,i]=stats(K60@v)[0]; Z["75"][j,i]=stats(K75@v)[0]
f,axs=plt.subplots(1,3,figsize=(17,5))
for ax,(k,t) in zip(axs,(("raw","raw signal"),("60","60° cone-averaged"),("75","75° cone-averaged"))):
    z=Z[k]; pc=ax.pcolormesh(X,Y,z,vmin=0,vmax=100,cmap="RdBu_r",shading="auto")
    if z.min()<50<z.max(): cs=ax.contour(X,Y,z,[50],colors="k",linewidths=1.5); ax.clabel(cs,fmt="50%%",fontsize=8)
    ax.set_xlabel("R / D  (galaxy distance / our distance to CoE)"); ax.set_ylabel("|v| / v_esc(D)")
    ax.set_title("%s\n(range %.1f–%.1f%% receding)"%(t,z.min(),z.max()),fontsize=10)
    print(k,"range %.1f–%.1f"%(z.min(),z.max()), "excluding y=0 row: %.1f–%.1f"%(z[1:].min(),z[1:].max()))
plt.colorbar(pc,ax=axs,label="% of sky receding (solid-angle weighted)")
f.savefig("cone_heat_recede_RD_vvesc.png",dpi=110,bbox_inches="tight")
