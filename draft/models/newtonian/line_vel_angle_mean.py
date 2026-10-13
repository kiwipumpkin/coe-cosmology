import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("ls_vel_angle.py").read().split("VB=")[0])
VB=0.4*vesc(M0,D0)
def plot(items,title,fn,cmap="viridis",cols=None):
    f,ax=plt.subplots(figsize=(9,5)); ax.axhline(0,color="k",lw=1.2)
    cs=cols or plt.get_cmap(cmap)(np.linspace(0.1,0.9,len(items)))
    W=np.abs(st)
    for (y,lab),c in zip(items,cs):
        ok=~np.isnan(y); m=(W[ok]*y[ok]).sum()/W[ok].sum(); inv=" *" if (~ok).any() else ""
        ax.plot(ang,y,color=c,lw=2,label=f"{lab}: mean {m:+,.0f}{inv}"); ax.axhline(m,color=c,lw=1.3,ls="--")
        print(fn,lab,round(m),"invalid" if inv else "")
    for x in (90,180,270): ax.axvline(x,color="0.8",lw=.7,ls=":")
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11)
    ax.set_ylabel("signed LOS relative velocity (km/s, + = receding)",fontsize=11); ax.set_title(title+"\n(dashed = solid-angle sky mean over valid directions; * = has invalid gaps)",fontsize=10.5)
    ax.legend(frameon=False,fontsize=9.5,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.13)); f.savefig(fn.replace(".png","_mean.png"),dpi=130,bbox_inches="tight"); plt.close(f)
xs=[-2,-1,-0.6,-0.4,0.4,1,2]
cols=[plt.cm.Blues(v) for v in (0.95,0.8,0.65,0.5)]+[plt.cm.Reds(v) for v in (0.5,0.7,0.9)]
plot([(vl(D0,R0,M0,x*vesc(M0,D0)),f"field {x:+} × v_esc ({x*vesc(M0,D0)/1e3:+,.0f} km/s)") for x in xs],
     "Vary velocity field (our shell speed / escape speed, all shells same energy)\nD = 58.2, R = 14.3 Gpc, M = 7.5×10⁵² kg; blue = collapsing, red = expanding","line_vel_angle_field.png",cols=cols)
plot([(vl(D0,R,M0,-VB),f"R = {R} Gpc") for R in (1,5,10,14.3,20,25)],"Vary galaxy distance R — collapsing field 0.4 × v_esc (−29,865 km/s), D = 58.2, M = 7.5×10⁵² kg","line_vel_angle_R.png")
plot([(vl(D,R0,M0,-VB),f"D = {D} Gpc") for D in (30,40,58.2,70,80)],"Vary our distance D — collapsing field −29,865 km/s fixed, R = 14.3, M = 7.5×10⁵² kg","line_vel_angle_D.png")
plot([(vl(D0,R0,M0*k,-VB),f"M = {k}× ({k*7.5:g}×10⁵² kg)") for k in (0.5,1,2,4)],"Vary mass M — collapsing field −29,865 km/s fixed, D = 58.2, R = 14.3","line_vel_angle_M.png")
