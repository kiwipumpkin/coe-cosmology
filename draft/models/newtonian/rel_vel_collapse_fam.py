import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_angle_cc.py").read().split("res={}")[0])
GM0=GM; D0=D; s0=s
def curve(gm=GM0,Dv=D0,sv=s0):
    xo=np.array([Dv,0.]); xp=xo+sv*n; rp=np.linalg.norm(xp,axis=1)
    vf=lambda r:-np.sqrt(2*gm/(r*GPC))
    vr=vf(rp)[:,None]*xp/rp[:,None]-vf(np.array(Dv))*xo/Dv; return (vr*n).sum(1)/1e3
w=np.abs(np.sin(t))
fams=[("M",[0.5,1,2,4,8],lambda x:curve(gm=GM0*x),lambda x:"M = %g× 7.5×10⁵² kg"%x,"mass"),
      ("D",[30,40,58.2,70,80],lambda x:curve(Dv=x),lambda x:"D = %g Gpc"%x,"observer distance D"),
      ("s",[1,5,10,14.3,20,25],lambda x:curve(sv=x),lambda x:"galaxy dist = %g Gpc"%x,"galaxy distance")]
base=curve()
for tag,vals,fn,lab,what in fams:
    f,ax=plt.subplots(figsize=(9,5.5)); ax.axhline(0,color="k",lw=1)
    cols=plt.cm.viridis(np.linspace(0,.9,len(vals)))
    for x,cl in zip(vals,cols):
        y=fn(x); m=(y*w).sum()/w.sum(); print(tag,x,round(m), "maxdiff_base" if np.allclose(y,base) else "")
        ax.plot(ang,y,color=cl,lw=2.2,label=lab(x)+"   (sky mean %+.0f km/s)"%m); ax.axhline(m,color=cl,ls="--",lw=1.3)
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE")
    ax.set_ylabel("line-of-sight relative velocity (km/s, + = receding)")
    ax.set_title("Collapsing shells (free fall): LOS velocity vs angle, varying "+what+"\nbaseline D = 58.2, s = 14.3 Gpc, M = 7.5×10⁵² kg; dashed = sky mean")
    ax.legend(frameon=False,ncol=2,fontsize=9,loc="upper center",bbox_to_anchor=(0.5,-0.13))
    f.savefig("rel_vel_collapse_fam_%s.png"%tag,dpi=130,bbox_inches="tight"); plt.close(f)
