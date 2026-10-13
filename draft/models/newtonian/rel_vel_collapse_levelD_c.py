import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_vel_collapse_combo_c.py").read().split("sets=")[0])
M,R,RST=2.0,10.0,60.0   # fixed absolute; D+R < RST so no line of sight crosses the stall radius (no kinks)
Ds=[20,25,30,35,40,45,49]
f,ax=plt.subplots(figsize=(9,6)); ax.axhline(0,color="k",lw=1); cols=plt.cm.viridis(np.linspace(0,.92,len(Ds)))
for D,cl in zip(Ds,cols):
    y=curve(M,D,R,RST); m=(y*w).sum()/w.sum(); fr=100*(w*(y>0)).sum()/w.sum()
    vus=np.sqrt(2*GM0*M/GPC*(1/D-1/RST))/c; ve=np.sqrt(2*GM0*M/(D*GPC))/c
    lab="D = %g Gpc | infall %.3f c (%.2f v_esc) | mean %+.3g c | rec %.0f%%"%(D,vus,vus/ve,m,fr); print(lab)
    ax.plot(ang,y,color=cl,lw=2,label=lab)
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE")
ax.set_ylabel("v_LOS / c   (+ = receding)")
ax.set_title("Collapsing shells (free fall from r_stall = 60 Gpc), level sets of observer distance D\nM = 2× 7.5×10⁵² kg, galaxy distance R = 10 Gpc (fixed); D + R < r_stall so no kinks",fontsize=10)
ax.legend(frameon=False,fontsize=8.5,loc="upper center",bbox_to_anchor=(0.5,-0.12))
f.savefig("rel_vel_collapse_levelD_c.png",dpi=130,bbox_inches="tight")
