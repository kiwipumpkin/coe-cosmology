import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_vel_collapse_combo_c.py").read().split("sets=")[0])
M,R,RST=2.0,5.0,56.0   # absolute; D+R < r_stall for all D -> no line of sight crosses r_stall (no kinks)
Ds=[20,30,40,45,48,50]; cols=plt.cm.viridis(np.linspace(0,.92,len(Ds)))
res={}
for sg,nm in ((-1,"collapse"),(1,"expand")):
    rows=[]
    for D in Ds:
        y=curve(M,D,R,RST) if sg==-1 else -curve(M,D,R,RST)   # expanding = mirror (outward, same |v| law)
        m=(y*w).sum()/w.sum(); fr=100*(w*(y>0)).sum()/w.sum(); vus=np.sqrt(2*GM0*M/GPC*(1/D-1/RST))/c; ve=np.sqrt(2*GM0*M/(D*GPC))/c
        rows.append((D,y,vus,vus/ve,m,fr)); print(nm,"D=%g v=%.3fc (%.2f v_esc) mean %+.4f rec %.1f%%"%(D,vus,vus/ve,m,fr))
    res[nm]=rows
lo=min(r[1].min() for v in res.values() for r in v); hi=max(r[1].max() for v in res.values() for r in v); pad=.05*(hi-lo)
ttl={"collapse":"COLLAPSING (free fall inward from r_stall = 56 Gpc)","expand":"EXPANDING (mirror: same speeds, outward)"}
def draw(ax,nm):
    ax.axhline(0,color="k",lw=1)
    for (D,y,vus,fv,m,fr),cl in zip(res[nm],cols):
        ax.plot(ang,y,color=cl,lw=2,label="D = %g Gpc | %s %.3f c (%.2f v_esc) | mean %+.3g c | rec %.0f%%"%(D,"infall" if nm=="collapse" else "outflow",vus,fv,m,fr))
    ax.set_ylim(lo-pad,hi+pad); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE")
    ax.set_ylabel("v_LOS / c   (+ = receding)"); ax.set_title(ttl[nm]+": level sets of observer distance D\nM = 2× 7.5×10⁵² kg, R = 5 Gpc fixed; D + R < r_stall so no kinks",fontsize=10)
    ax.legend(frameon=False,fontsize=8.3,loc="upper center",bbox_to_anchor=(0.5,-0.12))
for nm,fn in (("collapse","rel_vel_collapse_levelD_c.png"),("expand","rel_vel_expand_levelD_c.png")):
    f,ax=plt.subplots(figsize=(9,6)); draw(ax,nm); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
f,axs=plt.subplots(1,2,figsize=(18,6)); draw(axs[0],"expand"); draw(axs[1],"collapse"); f.savefig("rel_vel_levelD_both_c.png",dpi=110,bbox_inches="tight")
