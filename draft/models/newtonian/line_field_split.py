import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("ls_vel_angle.py").read().split("VB=")[0])
W=np.abs(st); ve=vesc(M0,D0)
def plot(xs,cmap,title,fn):
    f,ax=plt.subplots(figsize=(9,5)); ax.axhline(0,color="k",lw=1.2)
    cs=plt.get_cmap(cmap)(np.linspace(0.45,0.95,len(xs))); inv=np.zeros_like(ang,bool)
    for x,c in zip(xs,cs):
        y=vl(D0,R0,M0,x*ve); ok=~np.isnan(y); inv|=~ok; m=(W[ok]*y[ok]).sum()/W[ok].sum()
        ax.plot(ang,y,color=c,lw=2.2,label=f"field {x:+} × v_esc ({x*ve/1e3:+,.0f} km/s): mean {m:+,.0f}{' *' if (~ok).any() else ''}")
        ax.axhline(m,color=c,lw=1.3,ls="--"); print(fn,x,round(m))
    if inv.any():
        for a,b in ((0,ang[inv&(ang<180)].max()),(ang[inv&(ang>180)].min(),360)): ax.axvspan(a,b,color="0.9",zorder=0)
        ax.text(0.5,0.02,"grey span: some curves invalid (v² < 0, shells can't share our energy); * mean over valid directions only",transform=ax.transAxes,ha="center",fontsize=8.5,color="0.35")
    for x in (90,180,270): ax.axvline(x,color="0.8",lw=.7,ls=":")
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11)
    ax.set_ylabel("signed LOS relative velocity (km/s, + = receding)",fontsize=11)
    ax.set_title(title+"\n(dashed = solid-angle sky mean over valid directions; * = has invalid gaps)",fontsize=10.5)
    ax.legend(frameon=False,fontsize=9.5,ncol=1,loc="upper center",bbox_to_anchor=(0.5,-0.13)); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
plot([-0.4,-0.6,-1,-2],"Blues","Collapsing velocity fields (all shells same energy), D = 58.2, R = 14.3 Gpc, M = 7.5×10⁵² kg","line_vel_angle_field_col_mean.png")
plot([0.4,0.6,1,2],"Reds","Expanding velocity fields (all shells same energy), D = 58.2, R = 14.3 Gpc, M = 7.5×10⁵² kg","line_vel_angle_field_exp_mean.png")
