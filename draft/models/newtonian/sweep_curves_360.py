import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import importlib.util,sys
src=open("sweep_vel.py").read(); src=src[:src.index("S=np.linspace")]; exec(src)
th=np.radians(np.linspace(0,360,1441)); ct,st=np.cos(th),np.sin(th)
def vl(vf,D,s):
    x=D+s*ct; y=s*st; r=np.hypot(x,y); return vf(r)*(x*ct+y*st)/r-vf(np.array(D))*ct
ang=np.degrees(th); f,ax=plt.subplots(figsize=(9,4.8)); ax.axhline(0,color="k",lw=1.2)
for vf,col,lab,pc in ((vff(72.5),"#c0392b","free fall, r_max = 72.5 Gpc (stall at outer shell)",48),(vff(np.inf),"#2e5fa8","free fall, r_max = ∞",18),(vexp(1),"#7f3b3b","expanding, run-2 law, 1× mass",100)):
    v=vl(vf,58.2,14.3); ax.plot(ang,v/np.abs(v).max(),color=col,lw=2.2,label=f"{lab}, s = 14.3  (receding {pc}% of sky)")
lo,hi=ax.get_ylim(); ax.axhspan(0,hi,color="#fbe3e3",zorder=0); ax.axhspan(lo,0,color="#e3ebfb",zorder=0); ax.set_ylim(lo,hi)
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11)
ax.set_ylabel("LOS relative velocity / max |v|  (+ = receding)",fontsize=11); ax.set_title("Velocity vs angle at the extremes, D = 58.2 Gpc",fontsize=12)
ax.legend(frameon=False,fontsize=9.5,loc="upper center",bbox_to_anchor=(0.5,-0.14)); f.savefig("sweep_vel_curves_360.png",dpi=130,bbox_inches="tight")
