import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_angle_cc.py").read().split("res={}")[0])
def vlos(vf):
    vr=vf(rp)[:,None]*xp/rp[:,None]-vf(np.array(D))*xo/D; return (vr*n).sum(1)/1e3
ff=vlos(lambda r:-np.sqrt(2*GM/(r*GPC))); nm=vlos(lambda r:-k*r)
w=np.abs(np.sin(t)); mean=(ff*w).sum()/w.sum()
f,ax=plt.subplots(figsize=(9,5)); ax.axhline(0,color="k",lw=1)
ax.plot(ang,nm,color="0.5",ls="--",lw=1.8,label="collapsing, no mass (V = −k r)")
ax.plot(ang,ff,color="k",lw=2.4,label="collapsing, point mass (free fall from rest)")
ax.fill_between(ang,ff,0,where=ff>0,color="#c0392b",alpha=.35,label="receding")
ax.fill_between(ang,ff,0,where=ff<0,color="#2e5fa8",alpha=.35,label="approaching")
ax.axhline(mean,color="#7b3fa0",ls="--",lw=1.8,label="sky mean (solid-angle weighted) = %+.0f km/s"%mean)
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE")
ax.set_ylabel("line-of-sight relative velocity (km/s, + = receding)")
ax.set_title("Collapsing shells: relative velocity vs angle, D = 58.2, s = 14.3 Gpc, M = 7.5×10⁵² kg")
ax.legend(frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.13)); f.savefig("rel_vel_collapse_only.png",dpi=130,bbox_inches="tight")
print(mean, ff[0],ff[1440],ff.min(),nm[0])
