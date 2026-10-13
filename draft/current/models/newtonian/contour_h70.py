import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("integ_h70.py").read().split("Ms=[0,1,4,10]")[0].replace("print(","0 and print("))
G=6.674e-11; M1kg=Mscale*7.5e52
def fields(ms,ks):  # vectorized RK4, returns rf,v arrays (n,3000) in Gpc, m/s
    ri=np.linspace(0.02,5,3000)*GPC; GM=(np.array(ms)*Mscale*GM1)[:,None]; kk=np.array(ks)[:,None]
    r=np.tile(ri,(len(ms),1)); v=kk*ri; N=4000; h=T/N; f=lambda r:-GM/np.maximum(r,1e18)**2
    for _ in range(N):
        a1,b1_=v,f(r); a2,b2=v+.5*h*b1_,f(r+.5*h*a1); a3,b3=v+.5*h*b2,f(r+.5*h*a2); a4,b4=v+h*b3,f(r+h*a3)
        r=r+h/6*(a1+2*a2+2*a3+a4); v=v+h/6*(b1_+2*b2+2*b3+b4)
    return r/GPC,v
def stats(rf,v,Rv,Dd=D):
    xo=np.array([Dd,0.]); F=lambda x:np.interp(x,rf,v); hh=0.002*Dd/D
    def vl(R):
        xp=xo+R*n; rp=np.linalg.norm(xp,axis=1); return ((F(rp)[:,None]*xp/rp[:,None]-F(Dd)*xo/Dd)*n).sum(1)/c
    y=(vl(Rv+hh)-vl(Rv-hh))/(2*hh)*KMSMPC; mn=(w*y).sum()/w.sum(); sd=np.sqrt((w*(y-mn)**2).sum()/w.sum())
    return mn,sd,y[720],y[360]
def valid(rf,v,Dd,Rv):
    sel=(rf>=Dd-Rv)&(rf<=Dd+Rv); return np.all(np.diff(rf)>0) and v[sel].min()>0 and v[sel].max()<c
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"xtick.labelsize":16,"ytick.labelsize":16})
names=[("mean","sky-mean H (km/s/Mpc)"),("std","sky std dev of H (km/s/Mpc)"),("180","H toward CoE, 180° (km/s/Mpc)"),("90","H at 90° (km/s/Mpc)")]
def draw(X,Y,Z,bad,xl,yl,tag,dots,txt):
    for i,(nm,tt) in enumerate(names):
        Zi=np.ma.array(Z[...,i],mask=bad); f,ax=plt.subplots(figsize=(8,6)); ax.set_facecolor("0.75")
        cf=ax.contourf(X,Y,Zi,levels=14,cmap="viridis"); cs=ax.contour(X,Y,Zi,levels=cf.levels[::2],colors="k",linewidths=1.5)
        ax.clabel(cs,fmt="%.0f",fontsize=14); f.colorbar(cf,ax=ax)
        for (x,y),l in zip(dots,txt): ax.plot(x,y,"o",color="w",mec="k",ms=11); ax.annotate(l,(x,y),xytext=(6,6),textcoords="offset points",fontsize=16,color="w",weight="bold")
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(tt); f.tight_layout(); f.savefig("contour_H_%s_%s.png"%(tag,nm),dpi=110); plt.close(f)
# Map A: x=R/D, y=v/v_esc at observer (vary mass, k fixed)
ms=np.concatenate([[0],np.geomspace(0.1,60,45)]); rf,V=fields(ms,[k]*len(ms))
vesc=lambda m:np.sqrt(2*G*m*M1kg/(D*GPC)) if m>0 else np.inf
q=np.array([np.interp(D,rf[i],V[i])/vesc(m) for i,m in enumerate(ms)]); ratio=1/q  # v_esc/v ... use v/v_esc on y
RD=np.linspace(0.02,0.5,40); Z=np.zeros((len(ms),len(RD),4)); bad=np.zeros(Z.shape[:2],bool)
for i in range(len(ms)):
    for j,x in enumerate(RD):
        bad[i,j]=not valid(rf[i],V[i],D,x*D); Z[i,j]=stats(rf[i],V[i],x*D)
# y axis: v_esc/v (0 = no mass, grows with mass) -> more readable
yA=1/np.where(np.isinf(q),np.inf,q); yA[0]=0
X,Y=np.meshgrid(RD,yA); Zs=Z.copy()
dots=[(0.246,np.interp(m,ms,yA)) for m in (1,4,10)]
draw(X,Y,Z,bad,"R / D  (galaxy distance / CoE distance)","v_esc / v  at observer","RD_vesc",dots,["1×","4×","10×"])
for m in (1,4,10): print("m=%g: v/vesc=%.3f vesc/v=%.3f"%(m,1/np.interp(m,ms,yA),np.interp(m,ms,yA)))
for i in range(0,len(ms),4): print("m=%.2f vesc/v=%.3f first bad R/D=%s mean@.246=%.1f sd=%.1f"%(ms[i],yA[i],RD[bad[i]].min() if bad[i].any() else None,*stats(rf[i],V[i],0.246*D)[:2]))
np.savez("contour_A.npz",RD=RD,ms=ms,yA=yA,Z=Z,bad=bad)
# Map B: x=mass (x1x), y=k/k_70 at R/D=0.246
mB=np.linspace(0,30,31); kB=np.linspace(0.3,3,28); ZB=np.zeros((len(kB),len(mB),4)); bB=np.zeros((len(kB),len(mB)),bool)
for a,kf in enumerate(kB):
    rf,V=fields(mB,[kf*k]*len(mB))
    for b in range(len(mB)): bB[a,b]=not valid(rf[b],V[b],D,RMAX); ZB[a,b]=stats(rf[b],V[b],RMAX)
X,Y=np.meshgrid(mB,kB); draw(X,Y,ZB,bB,"mass M (units of 1× run)","initial expansion k / k_70","Mk",[(1,1),(4,1),(10,1)],["1×","4×","10×"])
np.savez("contour_B.npz",mB=mB,kB=kB,Z=ZB,bad=bB)
