"""Relative LOS velocity and acceleration vs angle, Chris convention (0 deg = away from CoE, 180 = toward).
D=58.2, R=s=14.3 Gpc. Point mass M=7.52e52 kg (paper run-2 central mass, GM = mu c^2 d1).
 expanding+mass: paper run-2 law v^2 = b1^2 r0^2 + 2mu(1/r-1/r0) (r=3 r0, units c, d1=D/3)
 collapsing+mass: free fall from rest at infinity, v = -sqrt(2GM/r)
 no mass: V = +/- k r, k = v_obs(run-2)/D, no gravity
Rel. acceleration = d2s/dt2 = a_rel.n + (|v_rel|^2 - v_los^2)/s."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
G=6.674e-11; c=2.998e8; GPC=3.0857e25; D=58.2; s=14.3
b1,mu=0.5866,0.09332; d1=D/3; GM=mu*c**2*d1*GPC
def v2(r): r0=r/3/d1; x=r/d1; return c*np.sqrt(b1**2*r0**2+2*mu*(1/x-1/r0))
k=v2(D)/D   # m/s per Gpc
ang=np.linspace(0,360,2881); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1)   # +x = away from CoE
xo=np.array([D,0.]); xp=xo+s*n; rp=np.linalg.norm(xp,axis=1)
cases={"expanding, no mass (V = +k r)":(lambda r:k*r,0,"#e8a0a0","--"),
       "expanding, point mass (paper run-2 law)":(v2,GM,"#c0392b","-"),
       "collapsing, no mass (V = −k r)":(lambda r:-k*r,0,"#9fb6dd","--"),
       "collapsing, point mass (free fall from rest)":(lambda r:-np.sqrt(2*GM/(r*GPC)),GM,"#2e5fa8","-")}
res={}
for name,(vf,gm,col,ls) in cases.items():
    vr=vf(rp)[:,None]*xp/rp[:,None]-vf(np.array(D))*xo/D
    a=(-gm*xp/rp[:,None]**3+gm*xo/D**3)/GPC**2
    vl=(vr*n).sum(1); acc=(a*n).sum(1)+((vr**2).sum(1)-vl**2)/(s*GPC)
    res[name]=(vl/1e3,acc,col,ls)
    def fl(y): i=np.where(np.diff(np.sign(y))!=0)[0]; return [round(float(ang[j]),1) for j in i if ang[j]<=180]
    print(name,"| v km/s: 0°=%.0f 90°=%.0f 180°=%.0f min %.0f@%.0f flips %s"%(vl[0]/1e3,vl[720]/1e3,vl[1440]/1e3,vl.min()/1e3,ang[vl.argmin()],fl(vl)),
          "| a m/s2: 0°=%.2e 90°=%.2e 180°=%.2e flips %s"%(acc[0],acc[720],acc[1440],fl(acc)))
for key,yl,fn in ((0,"line-of-sight relative velocity (km/s, + = receding)","rel_vel_angle_cc.png"),(1,"relative acceleration d²s/dt² (m/s², + = away)","rel_acc_angle_cc.png")):
    f,ax=plt.subplots(figsize=(9,5)); ax.axhline(0,color="k",lw=1)
    for name,v in res.items(): ax.plot(ang,v[key],color=v[2],ls=v[3],lw=2.2,label=name)
    lo,hi=ax.get_ylim(); ax.axhspan(0,hi,color="#fbe3e3",zorder=0); ax.axhspan(lo,0,color="#e3ebfb",zorder=0); ax.set_ylim(lo,hi)
    for x in (90,270): ax.axvline(x,color="0.75",lw=.8,ls=":")
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11); ax.set_ylabel(yl,fontsize=11)
    ax.set_title(("Relative velocity" if key==0 else "Relative acceleration")+" vs angle, D = 58.2, s = 14.3 Gpc, M = 7.5×10⁵² kg",fontsize=12)
    ax.legend(frameon=False,fontsize=10,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.13)); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
print("k =",k,"m/s/Gpc; v_obs exp",v2(D)/c,"c; v_obs ff",np.sqrt(2*GM/(D*GPC))/c,"c")
