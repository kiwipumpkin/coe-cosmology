"""Real integration replacing the sheet's one-step bolt-on.
Shells start at t0 with v = ±k r about CoE (k = 0.00805 c/Gpc, accepted baseline), r'' = -GM/r^2 (point mass at CoE), integrated
for T = 173.6 Gyr = the step implied by the sheet (strength GM dt/(D^2 k R) = 0.247 at M=1x, D=58.2, R=14.3).
At t0+T: observer = shell with r = 58.2 Gpc; galaxies at distance R = 14.3 Gpc in each direction (simultaneous, no light-travel time).
v_rel = v(r_g) r̂_g - v(D) r̂_o ; A/A_all = |v_rel| / |v_rel|(M=0) ; accel = a_rel·n + (|v_rel|^2 - v_los^2)/R  (rel_angle_cc.py formula)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
c=2.998e8; GPC=3.0857e25; GYR=3.156e16; b1,mu=0.5866,0.09332; GM1=mu*c**2*(58.2/3)*GPC
k=c*np.sqrt(b1**2+2*mu*(1/3-1))/58.2/GPC; T=0.247*(58.2*GPC)**2*k*14.3*GPC/GM1
D=58.2; R=14.3; ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1)
def field(sg,m):
    ri=np.linspace(5,400,6000)*GPC; GM=m*GM1
    sol=solve_ivp(lambda tt,y:np.concatenate([y[len(ri):],-GM/np.maximum(y[:len(ri)],1e20)**2]),(0,T),np.concatenate([ri,sg*k*ri]),rtol=1e-10,atol=1e5,method="DOP853")
    rf,vf=sol.y[:len(ri),-1]/GPC,sol.y[len(ri):,-1]
    ok=rf>0.5; rf,vf=rf[ok],vf[ok]
    assert np.all(np.diff(rf)>0), "shell crossing"
    assert rf.min()<D-R and rf.max()>D+R
    return lambda r:np.interp(r,rf,vf), GM
def calc(sg,m):
    vf,GM=field(sg,m); vr=vf(rp)[:,None]*xp/rp[:,None]-vf(np.array(D))*xo/D
    vl=(vr*n).sum(1); a=(-GM*xp/rp[:,None]**3+GM*xo/D**3)/GPC**2
    acc=(a*n).sum(1)+((vr**2).sum(1)-vl**2)/(R*GPC)
    return vl/c, np.linalg.norm(vr,axis=1), acc*GYR/c, vf(np.array(D))/c
Ms=[0,1,2,4,8]; res={}
for sg in (1,-1):
    for m in Ms: res[sg,m]=calc(sg,m)
# find recession mass for collapse
for m in (12,16,20,24,32):
    try:
        v=calc(-1,m)[0]
        if v.max()>0: res[-1,m]=calc(-1,m); res[1,m]=calc(1,m); Ms.append(m); break
    except AssertionError as e: print("M",m,e); break
mean=lambda y:(w*y).sum()/w.sum(); rec=lambda y:100*(w*(y>0)).sum()/w.sum()
out={}
for (sg,m),(vl,vabs,acc,vus) in res.items():
    A=vabs/res[sg,0][1]; out[sg,m]=(vl,A,acc)
    print("%s M=%g v_us %+.4fc | vLOS 0°%+.4f 90°%+.4f 180°%+.4f mean %+.4f rec %.1f%% | A 0°%.3f 90°%.3f 180°%.3f mean %.4f | acc mean %+.2e"%("exp" if sg>0 else "col",m,vus,vl[0],vl[360],vl[720],mean(vl),rec(vl),A[0],A[360],A[720],mean(A),mean(acc)))
cols=plt.cm.viridis(np.linspace(0.05,0.92,len(Ms)-1))
for idx,tag,yl in ((0,"vlos","v_LOS / c  (+ = receding)"),(1,"aaall","A / A_all = |v_rel| / |v_rel|(M=0)"),(2,"acc","relative LOS acceleration (c / Gyr)")):
    lo=min(out[x][idx].min() for x in out); hi=max(out[x][idx].max() for x in out); pad=.05*(hi-lo)
    for sg,sn in ((1,"expand"),(-1,"collapse")):
        f,ax=plt.subplots(figsize=(9,5.3)); ax.axhline(1 if idx==1 else 0,color="0.6",lw=.8)
        for j,m in enumerate(Ms):
            y=out[sg,m][idx]; lab="mean %+.4g"%mean(y)+(" c, rec %.0f%%"%rec(y) if idx==0 else (" c/Gyr" if idx==2 else ""))
            if m==0: ax.plot(ang,y,"k--",lw=1.6,label="no mass | "+lab)
            else: ax.plot(ang,y,color=cols[j-1],lw=2,label="M = %g× 7.5×10⁵² kg | "%m+lab)
        ax.set_ylim(lo-pad,hi+pad); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE"); ax.set_ylabel(yl)
        ax.set_title("%s — integrated shells (r'' = −GM/r², start v = %sk r, T = %.0f Gyr)\nobserver at D = 58.2, R = 14.3 Gpc at final time"%("EXPANDING" if sg>0 else "COLLAPSING","+" if sg>0 else "−",T/GYR),fontsize=10)
        ax.legend(frameon=False,fontsize=8.5,loc="upper center",bbox_to_anchor=(.5,-.13))
        f.savefig("integ_%s_%s_c.png"%(tag,sn),dpi=125,bbox_inches="tight"); plt.close(f)
print("T Gyr",T/GYR)
