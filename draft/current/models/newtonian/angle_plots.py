"""Line plots vs angle from CoE. (a) point mass 7.52e52 kg with paper run-2 velocities (consistent);
(b) exp L=6.07 Gpc, M(<D)=(4/3)pi D^3 * 2.9e-27 kg/m^3, V = k r with k = v_obs/D (homologous).
Curves divided by s (a_los grows ~ linearly with s)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
G=6.674e-11; c=2.998e8; GPC=3.0857e25; D=58.2
b1,mu,se=0.5866,0.09332,3.0; d1=D/se; GMc=mu*c**2*d1*GPC
MD=4/3*np.pi*(D*GPC)**3*2.9e-27; L=6.07
def m(x): return 1-np.exp(-x)*(1+x+x*x/2)
def vpm(r): r0=r/se/d1; x=r/d1; return c*np.sqrt(b1**2*r0**2+2*mu*(1/x-1/r0))
GMf={"pm":lambda r:GMc+0*r,"exp":lambda r:G*MD*m(r/L)/m(D/L)}
vf={"pm":vpm,"exp":lambda r:vpm(D)*r/D}
ang=np.linspace(0,180,3601); th=np.radians(180-ang)        # th from outward axis
n=np.stack([np.sin(th),0*th,np.cos(th)],1); xo=np.array([0,0,D])
def calc(p,s,sg):
    xg=xo+s*n; rg=np.linalg.norm(xg,axis=1)
    a=-(GMf[p](rg)/(rg*GPC)**3)[:,None]*xg*GPC+(GMf[p](D)/(D*GPC)**3)*xo*GPC
    alos=(a*n).sum(1)
    v=sg*(vf[p](rg)[:,None]*xg/rg[:,None]-vf[p](D)*xo/D); vl=(v*n).sum(1)
    sdd=alos+((v**2).sum(1)-vl**2)/(s*GPC)
    dspeed=(v*a).sum(1)/np.linalg.norm(v,axis=1)               # d|v_rel|/dt
    return alos,sdd,dspeed
def flips(y): i=np.where(np.diff(np.sign(y))!=0)[0]; return [round(float(ang[k]),1) for k in i]
T={"pm":"(a) point mass at CoE (7.5×10⁵² kg)","exp":"(b) concentrated mass, exp L = 6.07 Gpc"}
cols={1:"#d62728",7:"#9467bd",14:"#1f77b4"}
def plot(p,key,ylabel,title,fn):
    fig,ax=plt.subplots(figsize=(8,4.6)); rep={}
    for s in (1,7,14):
        if key=="alos": y=calc(p,s,1)[0]
        else: y=calc(p,s,1 if key=="exp_" else -1)[2]
        y=y/s; f=flips(y); rep[s]=f
        ax.plot(ang,y,color=cols[s],lw=2,label=f"s = {s} Gpc  (flips {', '.join(f'{a:.0f}°' for a in f)})")
        for a in f: ax.axvline(a,color=cols[s],ls=":",lw=1)
    ax.axhline(0,color="k",lw=0.8); ax.set_xlim(0,180); ax.set_xticks(range(0,181,30))
    lo,hi=ax.get_ylim(); ax.axhspan(0,hi,color="#fbe3e3",zorder=0); ax.axhspan(lo,0,color="#e3ebfb",zorder=0); ax.set_ylim(lo,hi)
    ax.set_xlabel("angle from CoE direction (deg)",fontsize=11); ax.set_ylabel(ylabel,fontsize=11)
    ax.set_title(title,fontsize=12); ax.legend(fontsize=10,frameon=False)
    fig.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(fig); print(fn,rep)
for p in ("pm","exp"):
    plot(p,"alos","a_rel · n / s   (m s⁻² per Gpc)",f"Relative acceleration seen by observer\n{T[p]}   (red = away, blue = toward)",f"angle_alos_{p}.png")
    plot(p,"exp_","d|v_rel|/dt / s   (m s⁻² per Gpc)",f"Change of recession speed, EXPANDING\n{T[p]}   (red = speeding up)",f"angle_speed_expanding_{p}.png")
    plot(p,"col_","d|v_rel|/dt / s   (m s⁻² per Gpc)",f"Change of approach speed, COLLAPSING\n{T[p]}   (red = speeding up)",f"angle_speed_collapsing_{p}.png")
