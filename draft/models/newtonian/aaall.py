"""Chris sheet col AF: A/A_all = |±kR n + dt (g_gal - g_obs)| / (kR); V = ±k r about CoE; point mass at CoE.
Strength GM dt/(D^2 k R) = 0.247 at D=58.2, R=14.3 (sheet value); GM dt/k held fixed when R or D vary. 0 deg = away from CoE."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D0,R0=58.2,14.3; C0=0.247*D0**2*R0
ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); W=np.abs(np.sin(t))
def A(sg=1,mf=1,R=R0,D=D0):
    xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1)
    dg=-xp/rp[:,None]**3+xo/D**3; return np.linalg.norm(sg*R*n+mf*C0*dg,axis=1)/R
def plot(items,title,fn,cols=None,cmap="viridis"):
    f,ax=plt.subplots(figsize=(9,5)); ax.axhline(1,color="0.4",lw=1.4,ls=":",label="no mass: 1 everywhere")
    cs=cols or plt.get_cmap(cmap)(np.linspace(0.1,0.9,len(items)))
    for (y,lab),c in zip(items,cs):
        m=(W*y).sum()/W.sum(); ax.plot(ang,y,color=c,lw=2,label=f"{lab}: sky mean {m:.3f}"); ax.axhline(m,color=c,lw=1.2,ls="--")
        print(fn,lab,"mean %.3f max %.3f@%.0f min %.3f@%.0f"%(m,y.max(),ang[y.argmax()],y.min(),ang[:721][y[:721].argmin()]))
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11)
    ax.set_ylabel("A / A_all = |±kR n + dt Δg| / kR",fontsize=11); ax.set_title(title+"\n(dashed = solid-angle sky mean)",fontsize=11)
    ax.legend(frameon=False,fontsize=9.5,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.13)); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
plot([(A(1),"expanding V = +k r"),(A(-1),"collapsing V = −k r")],"Chris sheet A/A_all: point mass at CoE (strength 0.247), D = 58.2, R = 14.3 Gpc","aaall_exp_vs_col.png",cols=["#c0392b","#2e5fa8"])
for sg,nm,w in ((1,"exp","expanding V = +k r"),(-1,"col","collapsing V = −k r")):
    plot([(A(sg,k),f"{k}× mass (strength {0.247*k:.3g})") for k in (0.5,1,2,4)],f"{w}: vary point mass, D = 58.2, R = 14.3","aaall_mass_"+nm+".png")
plot([(A(1,R=R),f"R = {R} Gpc (strength {C0/(D0**2*R):.3g})") for R in (1,5,10,14.3,20,25)],"expanding V = +k r: vary galaxy distance R, D = 58.2 (GM·dt/k fixed)","aaall_R.png")
plot([(A(1,D=D),f"D = {D} Gpc (strength {C0/(D**2*R0):.3g})") for D in (30,40,58.2,70,80)],"expanding V = +k r: vary our distance D, R = 14.3 (GM·dt/k fixed)","aaall_D.png")
