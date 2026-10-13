"""D level sets, fixed absolute params: M=2x7.5e52 kg, R=5 Gpc.
Collapse: free fall from rest at r_stall=56 Gpc, v=-sqrt(2GM(1/r-1/56)).
Expand: outward, energy-conserving, unbound: v=+sqrt(v_inf^2+2GM/r); v_inf anchored to accepted rel_angle_cc.py run-2
 expanding speed at D=58.2 (M=1 there): v_inf^2 = v2(58.2)^2 - 2GM1/58.2. (Run-2 law itself is imaginary for r<~52 Gpc at M=2.)
Rel. acceleration as rel_angle_cc.py: d2s/dt2 = a_rel.n + (|v_rel|^2 - v_los^2)/s. Unit c/Gyr."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
G=6.674e-11; c=2.998e8; GPC=3.0857e25; GYR=3.156e16
b1,mu=0.5866,0.09332; d1=58.2/3; GM1=mu*c**2*d1*GPC; GM=2*GM1
v2_582=c*np.sqrt(b1**2*(58.2/3/d1)**2+2*mu*(1/(58.2/d1)-1/(58.2/3/d1)))
vinf=np.sqrt(v2_582**2-2*GM1/(58.2*GPC)); RST=56.; R=5.
law={"collapse":lambda r:-np.sqrt(2*GM/GPC*np.maximum(0,1/r-1/RST)),"expand":lambda r:K_H*r}
ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
K_H=v2_582/58.2; C0=0.247*58.2**2*14.3   # accepted Chris-sheet/aaall.py: v_rel = kR n + dt*dg, strength GM dt/(D^2 kR)=0.247 at M=1,D=58.2,R=14.3, GM dt/k scales with M
def run(nm,D):
    if nm=="expand":
        xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1); dg=-xp/rp[:,None]**3+xo/D**3; dgn=(dg*n).sum(1)
        vl=K_H*R*(1+2*C0*dgn/R)                       # M=2x
        acc=GM*dgn/GPC**2                               # v_rel = kR n is radial -> transverse term 0
        return vl/c, acc*GYR/c
    vf=law[nm]; xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1)
    vr=vf(rp)[:,None]*xp/rp[:,None]-vf(np.array(D))*xo/D
    a=(-GM*xp/rp[:,None]**3+GM*xo/D**3)/GPC**2
    vl=(vr*n).sum(1); acc=(a*n).sum(1)+((vr**2).sum(1)-vl**2)/(R*GPC)
    return vl/c, acc*GYR/c
Lv=[20,30,40,45,48,50]; cols=plt.cm.viridis(np.linspace(0,.92,len(Lv)))
print("v_inf/c=%.4f"%(vinf/c))
R_={}
for nm in ("expand","collapse"):
    for D in Lv:
        v,a=run(nm,D); R_[nm,D]=(v,a)
        print("%-8s D=%g v_us=%+.3fc | mean v %+.4f rec %.1f%% | v@0 %+.4f @90 %+.4f @180 %+.4f | max %.4f@%.0f min %.4f@%.0f | acc mean %+.4f c/Gyr"%(nm,D,law[nm](D)/c,(w*v).sum()/w.sum(),100*(w*(v>0)).sum()/w.sum(),v[0],v[360],v[720],v.max(),ang[v.argmax()],v.min(),ang[v[:721].argmin()],(w*a).sum()/w.sum()))
def fam(k,ylab,fnames,tag):
    lo=min(R_[key][k].min() for key in R_); hi=max(R_[key][k].max() for key in R_); pad=.05*(hi-lo)
    axes={}
    def draw(ax,nm):
        ax.axhline(0,color="k",lw=1)
        for D,cl in zip(Lv,cols):
            y=R_[nm,D][k]; m=(w*y).sum()/w.sum(); fr=100*(w*(y>0)).sum()/w.sum()
            ax.plot(ang,y,color=cl,lw=2,label=("D = %g Gpc | v_us %+.3f c | mean %+.3g c | rec %.0f%%"%(D,law[nm](D)/c,m,fr)) if k==0 else ("D = %g Gpc | sky-mean accel %+.3g c/Gyr"%(D,m)))
        ax.set_ylim(lo-pad,hi+pad); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE"); ax.set_ylabel(ylab)
        ax.set_title(("EXPANDING (Hubble V = +k r, k = %.4f c/Gpc, decelerated by point mass; aaall.py signed LOS)"%(K_H/c) if nm=="expand" else "COLLAPSING (free fall from r_stall = 56 Gpc)")+" — %s, level sets of D\nM = 2× 7.5×10⁵² kg, R = 5 Gpc fixed"%tag,fontsize=10)
        ax.legend(frameon=False,fontsize=8.3,loc="upper center",bbox_to_anchor=(0.5,-0.12))
    for nm,fn in zip(("expand","collapse"),fnames[:2]):
        f,ax=plt.subplots(figsize=(9,6)); draw(ax,nm); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
    f,axs=plt.subplots(1,2,figsize=(18,6)); draw(axs[0],"expand"); draw(axs[1],"collapse"); f.savefig(fnames[2],dpi=110,bbox_inches="tight"); plt.close(f)
fam(0,"v_LOS / c   (+ = receding)",["rel_vel_expand_levelD_c.png","rel_vel_collapse_levelD_c.png","rel_vel_levelD_both_c.png"],"LOS velocity")
fam(1,"relative LOS acceleration d²s/dt²  (c / Gyr, + = away)",["rel_acc_levelD_expand_c.png","rel_acc_levelD_collapse_c.png","rel_acc_levelD_both_c.png"],"relative acceleration")
Ds=np.linspace(15,50.9,120); f,(a1,a2,a3)=plt.subplots(3,1,figsize=(8.5,9.5),sharex=True)
for nm,col,lab in (("expand","#c0392b","expanding (V = +k r with point-mass deceleration)"),("collapse","#2e5fa8","collapsing (free fall from 56 Gpc)")):
    S=np.array([[ (w*vv).sum()/w.sum(),100*(w*(vv>0)).sum()/w.sum(),(w*aa).sum()/w.sum()] for vv,aa in (run(nm,D) for D in Ds)])
    P=np.array([[ (w*R_[nm,D][0]).sum()/w.sum(),100*(w*(R_[nm,D][0]>0)).sum()/w.sum(),(w*R_[nm,D][1]).sum()/w.sum()] for D in Lv])
    for ax,j in ((a1,0),(a2,1),(a3,2)): ax.plot(Ds,S[:,j],color=col,lw=2.2,label=lab); ax.plot(Lv,P[:,j],"o",color=col,mec="k")
    if nm=="collapse":
        i=np.where(np.diff(np.sign(S[:,0])))[0]; j=np.where(np.diff(np.sign(S[:,1]-50)))[0]; print("collapse mean0 D≈",Ds[i] ,"50% D≈",Ds[j])
a1.axhline(0,color="k",lw=1); a2.axhline(50,color="k",ls=":"); a3.axhline(0,color="k",lw=1)
a1.set_ylabel("sky-average v_LOS / c"); a2.set_ylabel("% of sky receding"); a3.set_ylabel("sky-average rel. accel\n(c / Gyr)"); a3.set_xlabel("observer distance to CoE, D (Gpc)")
a1.legend(frameon=False,fontsize=9); a1.set_title("Sky averages vs D (solid-angle weighted; dots = level-set D)\nM = 2× 7.5×10⁵² kg, R = 5 Gpc",fontsize=10)
f.tight_layout(); f.savefig("skymean_vs_D_c.png",dpi=120)
f,ax=plt.subplots(figsize=(8.5,5)); ax.axhline(0,color="k",lw=1)
for nm,col,lab in (("expand","#c0392b","expanding (V = +k r with point-mass deceleration)"),("collapse","#2e5fa8","collapsing (free fall from 56 Gpc)")):
    A=[(w*run(nm,D)[1]).sum()/w.sum() for D in Ds]; ax.plot(Ds,A,color=col,lw=2.2,label=lab)
    ax.plot(Lv,[(w*R_[nm,D][1]).sum()/w.sum() for D in Lv],"o",color=col,mec="k")
ax.set_xlabel("observer distance to CoE, D (Gpc)"); ax.set_ylabel("sky-average relative LOS acceleration (c / Gyr)")
ax.set_title("All-sky (solid-angle weighted) relative acceleration vs D (dots = level-set D)\nM = 2× 7.5×10⁵² kg, R = 5 Gpc",fontsize=10); ax.legend(frameon=False)
f.savefig("skyacc_vs_D_c.png",dpi=130,bbox_inches="tight")
