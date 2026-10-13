"""D level sets with the accepted Chris-sheet method (aaall.py) for BOTH cases:
v_rel = ±kR n + dt*Δg ; signed LOS = n·v_rel = ±kR + dt Δg·n ; A/A_all = |v_rel|/kR.
M = 1× (sheet strength GM dt/(D²kR) = 0.247 at D=58.2, R=14.3), R = 14.3 Gpc, GM·dt/k fixed. k from rel_angle_cc.py (v_obs/D at 58.2).
Accel (rel_angle_cc.py formula): a_rel·n + (|v_rel|²-v_los²)/s with v_rel = ±kR n (radial) -> GM Δg·n. Unit c/Gyr."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
c=2.998e8; GPC=3.0857e25; GYR=3.156e16; b1,mu=0.5866,0.09332; d1=58.2/3; GM=mu*c**2*d1*GPC
k=c*np.sqrt(b1**2+2*mu*(1/3-1))/58.2          # m/s per Gpc  (v2(58.2)/58.2)
R=14.3; C0=0.247*58.2**2*14.3
ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
def run(sg,D):
    xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1); dg=-xp/rp[:,None]**3+xo/D**3
    vrel=sg*R*n+C0*dg                              # in units of k (Gpc)
    return k*(vrel*n).sum(1)/c, np.linalg.norm(vrel,axis=1)/R, GM*(dg*n).sum(1)/GPC**2*GYR/c
Lv=[40,45,50,58.2,70,80]; cols=plt.cm.viridis(np.linspace(0,.92,len(Lv))); res={(s,D):run(s,D) for s in (1,-1) for D in Lv}
nm={1:"EXPANDING (V = +k r)",-1:"COLLAPSING (V = −k r)"}; mean=lambda y:(w*y).sum()/w.sum(); rec=lambda y:100*(w*(y>0)).sum()/w.sum()
for s in (1,-1):
    for D in Lv:
        v,A,a=res[s,D]; print("%s D=%g | vLOS 0°%+.4f 90°%+.4f 180°%+.4f mean %+.4f rec %.0f%% | A 0°%.3f 90°%.3f 180°%.3f max %.3f@%.0f min %.3f@%.0f mean %.3f | acc mean %+.2e"%(nm[s][:4],D,v[0],v[360],v[720],mean(v),rec(v),A[0],A[360],A[720],A.max(),ang[A[:721].argmax()],A.min(),ang[A[:721].argmin()],mean(A),mean(a)))
def fam(idx,ylab,files,tag,legf):
    lo=min(res[x][idx].min() for x in res); hi=max(res[x][idx].max() for x in res); pad=.05*(hi-lo)
    def draw(ax,s):
        ax.axhline(1 if idx==1 else 0,color="k",lw=1)
        for D,cl in zip(Lv,cols): y=res[s,D][idx]; ax.plot(ang,y,color=cl,lw=2.4 if D==58.2 else 1.8,label=legf(D,y))
        ax.set_ylim(lo-pad,hi+pad); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE"); ax.set_ylabel(ylab)
        ax.set_title(nm[s]+" — "+tag+", level sets of D\npoint mass 1× 7.5×10⁵² kg (sheet strength 0.247 at D=58.2), R = 14.3 Gpc; thick = accepted baseline",fontsize=10)
        ax.legend(frameon=False,fontsize=8.5,loc="upper center",bbox_to_anchor=(0.5,-0.12),ncol=2)
    for s,fn in zip((1,-1),files[:2]):
        f,ax=plt.subplots(figsize=(9,6)); draw(ax,s); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
    f,axs=plt.subplots(1,2,figsize=(18,6)); draw(axs[0],1); draw(axs[1],-1); f.savefig(files[2],dpi=110,bbox_inches="tight"); plt.close(f)
fam(0,"v_LOS / c  (+ = receding)",["rel_vel_expand_levelD_v2_c.png","rel_vel_collapse_levelD_v2_c.png","rel_vel_levelD_both_v2_c.png"],"signed LOS velocity",lambda D,y:"D = %g | mean %+.3g c | rec %.0f%%"%(D,mean(y),rec(y)))
fam(1,"A / A_all = |±kR n + dt Δg| / kR",["aaall_levelD_expand.png","aaall_levelD_collapse.png","aaall_levelD_both.png"],"A/A_all",lambda D,y:"D = %g | sky mean %.3f"%(D,mean(y)))
fam(2,"relative LOS acceleration (c / Gyr, + = away)",["rel_acc_expand_levelD_v2_c.png","rel_acc_collapse_levelD_v2_c.png","rel_acc_levelD_both_v2_c.png"],"relative acceleration",lambda D,y:"D = %g | sky mean %+.2g c/Gyr"%(D,mean(y)))
Ds=np.linspace(38,85,150)
f,(a1,a2,a3)=plt.subplots(3,1,figsize=(8.5,9),sharex=True)
for s,col in ((1,"#c0392b"),(-1,"#2e5fa8")):
    S=np.array([[mean(r[0]),rec(r[0]),mean(r[1])] for r in (run(s,D) for D in Ds)]); P=np.array([[mean(res[s,D][0]),rec(res[s,D][0]),mean(res[s,D][1])] for D in Lv])
    for ax,j in ((a1,0),(a2,1),(a3,2)): ax.plot(Ds,S[:,j],color=col,lw=2.2,label=nm[s]); ax.plot(Lv,P[:,j],"o",color=col,mec="k")
a1.axhline(0,color="k"); a2.axhline(50,color="k",ls=":"); a3.axhline(1,color="k",ls=":")
a1.set_ylabel("sky-average v_LOS / c"); a2.set_ylabel("% sky receding"); a3.set_ylabel("sky-average A/A_all"); a3.set_xlabel("observer distance to CoE, D (Gpc)")
a1.legend(frameon=False); a1.set_title("Sky averages vs D (solid-angle weighted; dots = level sets)\nM = 1×, R = 14.3 Gpc, accepted sheet method",fontsize=10); f.tight_layout(); f.savefig("skymean_vs_D_v2_c.png",dpi=120); plt.close(f)
f,ax=plt.subplots(figsize=(8.5,5)); ax.axhline(0,color="k")
for s,col,ls in ((1,"#c0392b","-"),(-1,"#2e5fa8","--")):
    ax.plot(Ds,[mean(run(s,D)[2]) for D in Ds],color=col,ls=ls,lw=2.2,label=nm[s]); ax.plot(Lv,[mean(res[s,D][2]) for D in Lv],"o",color=col,mec="k")
ax.set_xlabel("observer distance to CoE, D (Gpc)"); ax.set_ylabel("sky-average relative LOS acceleration (c / Gyr)"); ax.legend(frameon=False)
ax.set_title("All-sky relative acceleration vs D (identical for both: gravity term only)\nM = 1×, R = 14.3 Gpc",fontsize=10); f.savefig("skyacc_vs_D_v2_c.png",dpi=130,bbox_inches="tight")
