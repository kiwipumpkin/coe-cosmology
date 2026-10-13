"""D level sets with the accepted Chris-sheet method (aaall.py) for BOTH cases:
v_rel = ±kR n + dt*Δg ; signed LOS = n·v_rel = ±kR + dt Δg·n ; A/A_all = |v_rel|/kR.
M = 1× (sheet strength GM dt/(D²kR) = 0.247 at D=58.2, R=14.3), R = 14.3 Gpc, GM·dt/k fixed. k from rel_angle_cc.py (v_obs/D at 58.2).
Accel (rel_angle_cc.py formula): a_rel·n + (|v_rel|²-v_los²)/s with v_rel = ±kR n (radial) -> GM Δg·n. Unit c/Gyr."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
c=2.998e8; GPC=3.0857e25; GYR=3.156e16; b1,mu=0.5866,0.09332; d1=58.2/3; GM=mu*c**2*d1*GPC
k=c*np.sqrt(b1**2+2*mu*(1/3-1))/58.2          # m/s per Gpc  (v2(58.2)/58.2)
R=14.3; C0=0.247*58.2**2*14.3
ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
def run(sg,D,mf=1.0):
    xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1); dg=-xp/rp[:,None]**3+xo/D**3
    vrel=sg*R*n+mf*C0*dg                              # in units of k (Gpc)
    return k*(vrel*n).sum(1)/c, np.linalg.norm(vrel,axis=1)/R, GM*(dg*n).sum(1)/GPC**2*GYR/c
Lv=[40,45,50,58.2,70,80]; cols=plt.cm.viridis(np.linspace(0,.92,len(Lv))); res={(s,D):run(s,D) for s in (1,-1) for D in Lv}
nm={1:"EXPANDING (V = +k r)",-1:"COLLAPSING (V = −k r)"}; mean=lambda y:(w*y).sum()/w.sum(); rec=lambda y:100*(w*(y>0)).sum()/w.sum()
for s in ():
    for D in Lv:
        v,A,a=res[s,D]; print("%s D=%g | vLOS 0°%+.4f 90°%+.4f 180°%+.4f mean %+.4f rec %.0f%% | A 0°%.3f 90°%.3f 180°%.3f max %.3f@%.0f min %.3f@%.0f mean %.3f | acc mean %+.2e"%(nm[s][:4],D,v[0],v[360],v[720],mean(v),rec(v),A[0],A[360],A[720],A.max(),ang[A[:721].argmax()],A.min(),ang[A[:721].argmin()],mean(A),mean(a)))
def fam(idx,ylab,files,tag,legf):
    lo=min(res[x][idx].min() for x in res); hi=max(res[x][idx].max() for x in res); pad=.05*(hi-lo)
    def draw(ax,s):
        ax.axhline(1 if idx==1 else 0,color="k",lw=1)
        for D,cl in zip(Lv,cols): y=res[s,D][idx]; ax.plot(ang,y,color=cl,lw=2.4 if D==58.2 else 1.8,label=legf(D,y))
        y0=run(s,58.2,0.0)[idx]; ax.plot(ang,y0,'k--',lw=1.5,label='no mass (%.4f)'%y0.mean())
        ax.set_ylim(lo-pad,hi+pad); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE"); ax.set_ylabel(ylab)
        ax.set_title(nm[s]+" — "+tag+", level sets of D\npoint mass 1× 7.5×10⁵² kg (sheet strength 0.247 at D=58.2), R = 14.3 Gpc; thick = accepted baseline",fontsize=10)
        ax.legend(frameon=False,fontsize=8.5,loc="upper center",bbox_to_anchor=(0.5,-0.12),ncol=2)
    for s,fn in zip((1,-1),files[:2]):
        pass
    f,axs=plt.subplots(1,2,figsize=(18,6)); draw(axs[0],1); draw(axs[1],-1); f.savefig(files[2],dpi=110,bbox_inches="tight"); plt.close(f)
fam(0,"v_LOS / c  (+ = receding)",["x","x","rel_vel_levelD_both_v2_c_nomass.png"],"signed LOS velocity",lambda D,y:"D = %g | mean %+.3g c | rec %.0f%%"%(D,mean(y),rec(y)))
fam(1,"A / A_all = |±kR n + dt Δg| / kR",["x","x","aaall_levelD_both_nomass.png"],"A/A_all",lambda D,y:"D = %g | sky mean %.3f"%(D,mean(y)))

for sg in (1,-1):
    for D in Lv:
        v,A,_=run(sg,D,0.0); print(sg,D,"v %.6f..%.6f A %.6f..%.6f"%(v.min(),v.max(),A.min(),A.max()))
