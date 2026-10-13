import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_vel_collapse_combo_c.py").read().split("sets=")[0])
def frac(m=1,Dv=58.2,sv=14.3,q=1.5):
    y=curve(m,Dv,sv,q*Dv); return 100*(w*(y>0)).sum()/w.sum()
Ms=np.geomspace(0.1,10,25); Qs=np.linspace(1.01,4,25); Ds=np.linspace(25,90,27); Rs=np.linspace(1,30,30)
maps=[("stallD_M",Ms,Qs,lambda a,b:frac(m=a,q=b),"central mass M (× 7.5×10⁵² kg)","stall radius / D",True),
      ("D_R",Rs,Ds,lambda a,b:frac(sv=a,Dv=b),"galaxy distance R (Gpc)","our distance to CoE D (Gpc)",False),
      ("M_D",Ds,Ms,lambda a,b:frac(Dv=a,m=b),"our distance to CoE D (Gpc)","central mass M (× 7.5×10⁵² kg)",False),
      ("R_stallD",Qs,Rs,lambda a,b:frac(q=a,sv=b),"stall radius / D","galaxy distance R (Gpc)",False)]
for tag,X,Y,fn,xl,yl,logx in maps:
    Z=np.array([[fn(a,b) for a in X] for b in Y]); print(tag,"%.1f–%.1f%%"%(Z.min(),Z.max()))
    f,ax=plt.subplots(figsize=(6.5,5)); pc=ax.pcolormesh(X,Y,Z,vmin=0,vmax=100,cmap="RdBu_r",shading="auto")
    if Z.min()<50<Z.max(): cs=ax.contour(X,Y,Z,[50],colors="k",linewidths=1.5); ax.clabel(cs,fmt="50%%")
    if logx: ax.set_xscale("log")
    if tag=="M_D": ax.set_yscale("log")
    ax.set_xlabel(xl); ax.set_ylabel(yl); plt.colorbar(pc,label="% of sky receding (solid-angle weighted)")
    ax.set_title("Collapse (free fall from stall radius): % sky receding\nothers at base: M=1×, D=58.2, R=14.3 Gpc, stall=1.5·D",fontsize=10)
    f.savefig("heat_recede_%s.png"%tag,dpi=120,bbox_inches="tight"); plt.close(f)
