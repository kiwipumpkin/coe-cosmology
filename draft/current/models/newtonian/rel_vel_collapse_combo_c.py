import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_angle_cc.py").read().split("res={}")[0])
GM0=GM; w=np.abs(np.sin(t))
def curve(m=1,Dv=58.2,sv=14.3,rst=np.inf):
    gm=GM0*m; xo=np.array([Dv,0.]); xp=xo+sv*n; rp=np.linalg.norm(xp,axis=1)
    vf=lambda r:-np.sqrt(2*gm/GPC*np.maximum(0,1/r-1/rst))   # free fall from rest at stall radius; shells beyond rst stalled (v=0)
    vr=vf(rp)[:,None]*xp/rp[:,None]-vf(np.array(Dv))*xo/Dv; return (vr*n).sum(1)/c
sets=[dict(m=1,Dv=58.2,sv=14.3,rst=np.inf),dict(m=1,Dv=58.2,sv=14.3,rst=58.2*1.02),dict(m=2,Dv=40,sv=8,rst=41),
      dict(m=1,Dv=58.2,sv=14.3,rst=58.2*1.5),dict(m=4,Dv=70,sv=25,rst=72),dict(m=0.5,Dv=30,sv=25,rst=np.inf),
      dict(m=4,Dv=40,sv=2,rst=40.05),dict(m=1,Dv=40,sv=20,rst=60)]
f,ax=plt.subplots(figsize=(9,6)); ax.axhline(0,color="k",lw=1)
cols=plt.cm.turbo(np.linspace(.05,.95,len(sets)))
for p,cl in zip(sets,cols):
    y=curve(**p); m=(y*w).sum()/w.sum(); fr=(w*(y>0)).sum()/w.sum()
    rs="∞" if np.isinf(p["rst"]) else "%g"%p["rst"]
    lab="M=%g×, D=%g, R=%g, r_stall=%s Gpc | mean %+.3g c, rec %.0f%%"%(p["m"],p["Dv"],p["sv"],rs,m,100*fr)
    print(lab, "v_ours/c=%.3g"%(np.sqrt(2*GM0*p["m"]/GPC*max(0,1/p["Dv"]-1/p["rst"]))/c))
    ax.plot(ang,y,color=cl,lw=2.6 if np.isinf(p["rst"]) and p["Dv"]==58.2 else 1.8,label=lab+(" (baseline)" if np.isinf(p["rst"]) and p["Dv"]==58.2 else ""))
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE")
ax.set_ylabel("v_LOS / c   (+ = receding)")
ax.set_title("Collapsing shells (free fall from stall radius), point mass at CoE: combined parameter changes\nM in units of 7.5×10⁵² kg; shells beyond r_stall are at rest")
ax.legend(frameon=False,ncol=1,fontsize=8.5,loc="upper center",bbox_to_anchor=(0.5,-0.12))
f.savefig("rel_vel_collapse_combo_c.png",dpi=130,bbox_inches="tight")
