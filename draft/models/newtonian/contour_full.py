import numpy as np
src=open("contour_h70.py").read(); exec(src.split("# Map A")[0])
_stats=stats
def stats(rf,v,Rv,Dd=D):
    o=np.argsort(rf); return _stats(rf[o],v[o],Rv,Dd)   # robust if shells cross
def flags(rf,v,Rv):
    sel=(rf>=D-Rv)&(rf<=D+Rv); return v[sel].max()>=c, (v[sel].min()<=0) or not np.all(np.diff(rf)>0)
def draw2(X,Y,Z,FC,FT,xl,yl,tag,dots,ylim=None):
    for i,(nm,tt) in enumerate(names):
        f,ax=plt.subplots(figsize=(8,6))
        cf=ax.contourf(X,Y,Z[...,i],levels=14,cmap="viridis"); cs=ax.contour(X,Y,Z[...,i],levels=cf.levels[::2],colors="k",linewidths=1.5)
        ax.clabel(cs,fmt="%.0f",fontsize=14); f.colorbar(cf,ax=ax)
        for F,lab in ((FC,"v > c"),(FT,"inner shells turn around")):
            if F.any() and not F.all():
                b=ax.contour(X,Y,F.astype(float),levels=[0.5],colors="w",linewidths=2,linestyles="--")
                ax.clabel(b,fmt={0.5:lab},fontsize=14,colors="w")
        for (x,y),l in zip(dots,["1×","4×","10×"]): ax.plot(x,y,"o",color="w",mec="k",ms=11); ax.annotate(l,(x,y),xytext=(6,6),textcoords="offset points",fontsize=16,color="w",weight="bold")
        if ylim: ax.set_ylim(*ylim)
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(tt); f.tight_layout(); f.savefig("contour_H_%s_%s.png"%(tag,nm),dpi=110); plt.close(f)
# Map A
ms=np.concatenate([[0],np.geomspace(0.1,15,60)]); rf,V=fields(ms,[k]*len(ms))
yA=np.array([0]+[np.sqrt(2*G*m*M1kg/(D*GPC))/np.interp(D,np.sort(rf[i]),V[i][np.argsort(rf[i])]) for i,m in enumerate(ms) if m>0])
keep=(yA>=0)&(yA<=4.2); print("Map A max mass",ms[keep].max())
RD=np.linspace(0.02,0.5,40); Z=np.zeros((len(ms),len(RD),4)); FC=np.zeros(Z.shape[:2],bool); FT=FC.copy()
for i in np.where(keep)[0]:
    for j,x in enumerate(RD): FC[i,j],FT[i,j]=flags(rf[i],V[i],x*D); Z[i,j]=stats(rf[i],V[i],x*D)
ms,yA,Z,FC,FT=ms[keep],yA[keep],Z[keep],FC[keep],FT[keep]
X,Y=np.meshgrid(RD,yA); draw2(X,Y,Z,FC,FT,"R / D  (galaxy distance / CoE distance)","v_esc / v  at observer","RD_vesc",[(0.246,np.interp(m,ms,yA)) for m in (1,4,10)],(0,yA.max()))
print("A ranges: mean %.0f..%.0f sd %.0f..%.0f"%(Z[...,0].min(),Z[...,0].max(),Z[...,1].min(),Z[...,1].max()))
import sys; sys.exit() if 0 else None
# Map B
mB=np.linspace(0,30,31); kB=np.linspace(0.3,3,28); ZB=np.zeros((len(kB),len(mB),4)); BC=np.zeros((len(kB),len(mB)),bool); BT=BC.copy()
for a,kf in enumerate(kB):
    rf,V=fields(mB,[kf*k]*len(mB))
    for b in range(len(mB)): BC[a,b],BT[a,b]=flags(rf[b],V[b],RMAX); ZB[a,b]=stats(rf[b],V[b],RMAX)
X,Y=np.meshgrid(mB,kB); draw2(X,Y,ZB,BC,BT,"mass M (units of 1× run)","initial expansion k / k_70","Mk",[(1,1),(4,1),(10,1)])
print("B ranges: mean %.0f..%.0f sd %.0f..%.0f  nan:%d"%(ZB[...,0].min(),ZB[...,0].max(),ZB[...,1].min(),ZB[...,1].max(),np.isnan(ZB).sum()))
