"""Free-fall collapse under point mass at CoE: v(r) = -sqrt(2GM(1/r-1/rmax)) radial. Observer D=58.2, galaxy at s=R=14.3 Gpc.
Mass scale: 1x = sheet strength 0.247 = GM dt/(D^2 k R), k=2.18e-18 /s, dt=1 s (sheet values) -> GM in SI.
Quantity: v_rel + dt*(g_gal-g_obs); A/A_all = |.|/sky-mean|.| ; signed = (.).n / sky-mean|.|."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
GPC=3.0857e25; c=2.998e8; D,R=58.2*GPC,14.3*GPC; k=2.18e-18; dt=1.0
GM1=0.247*D**2*k*R/dt
ang=np.linspace(0,360,2881); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1)
xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1)
w=np.abs(np.sin(t)); w/=w.sum()                     # sky (solid-angle) weight, axisymmetric about CoE axis
def field(GM,law,rmax=np.inf,grav=True):
    if law=="ff": vm=lambda r:-np.sqrt(2*GM*(1/r-1/rmax))
    else: vm=lambda r:k*r
    vr=vm(rp)[:,None]*xp/rp[:,None]-vm(D)*xo/D
    if grav: vr=vr+dt*(-GM*xp/rp[:,None]**3+GM*xo/D**3)
    mag=np.linalg.norm(vr,axis=1); m=(w*mag).sum(); return mag/m,(vr*n).sum(1)/m,mag,vm(D)
def desc(s):
    i=np.where(np.diff(np.sign(s))!=0)[0]; return dict(rec_frac=round(float((w*(s>0)).sum()),3),flips=[round(float(ang[j]),1) for j in i if ang[j]<=180],
        sign_0=round(float(s[0]),3),sign_90=round(float(s[720]),3),sign_180=round(float(s[1440]),3),min=round(float(s.min()),3),min_at=float(ang[s.argmin()]))
F=[0.5,1,2,4,8]; cols=plt.cm.viridis(np.linspace(0.2,0.9,len(F)))
rout=72.5*GPC
for rm,lab in ((np.inf,"rmax=inf"),(1.5*rout,"rmax=1.5x72.5"),(2*rout,"rmax=2x72.5")):
    for fm in F:
        A,S,mag,vo=field(GM1*fm,"ff",rm); _,S0,_,_=field(GM1*fm,"ff",rm,grav=False)
        print(lab,fm,"v_obs/c=%.3g"%(abs(vo)/c),"grav step / |v_rel| =%.1e"%(dt*GM1*fm/D**2/mag.mean()),desc(S),"Amax %.3f@%g Amin %.3f@%g"%(A.max(),ang[A.argmax()],A.min(),ang[A.argmin()]),"no-step same:",np.allclose(S,S0,atol=1e-6))
def plot(kind,fn):
    f,ax=plt.subplots(figsize=(9,5))
    ax.axhline(0 if kind=="S" else 1,color="0.5",lw=1.4,ls="--",label="0× : no motion (flat)" if kind=="A" else "zero")
    for fm,cl in zip(F,cols):
        A,S,_,_=field(GM1*fm,"ff"); ax.plot(ang,A if kind=="A" else S,color=cl,lw=5-0.7*F.index(fm),alpha=.9,label=f"{fm}× sheet mass, r_max = ∞")
    for rm,ls,l in ((1.5*rout,":","r_max = 1.5×72.5 Gpc (any mass)"),(2*rout,"-.","r_max = 2×72.5 Gpc (any mass)")):
        A,S,_,_=field(GM1,"ff",rm); ax.plot(ang,A if kind=="A" else S,color="k",lw=1.6,ls=ls,label=l)
    for x in (90,270): ax.axvline(x,color="0.75",lw=.8,ls=":")
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11)
    ax.set_ylabel("A / A_all  (|v_rel + dt Δg| / sky mean)" if kind=="A" else "signed LOS / sky-mean speed  (+ = receding)",fontsize=11)
    ax.set_title("Free-fall collapse, v = −√(2GM(1/r − 1/r_max)), point mass at CoE"+(", signed" if kind=="S" else ""),fontsize=12)
    ax.legend(frameon=False,fontsize=9,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.13)); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
plot("A","ff_col_masses_cc.png"); plot("S","ff_col_masses_signed_cc.png")
f,ax=plt.subplots(figsize=(9,5)); ax.axhline(0,color="k",lw=1)
_,Se,_,_=field(GM1,"exp"); _,Sc,_,_=field(GM1,"ff")
ax.plot(ang,Se,color="#c0392b",lw=2.2,label="expanding, Sedov–Taylor V = +k r, 1× mass"); ax.plot(ang,Sc,color="#2e5fa8",lw=2.2,label="free-fall collapse, r_max = ∞, 1× mass")
lo,hi=ax.get_ylim(); ax.axhspan(0,hi,color="#fbe3e3",zorder=0); ax.axhspan(lo,0,color="#e3ebfb",zorder=0); ax.set_ylim(lo,hi)
for x in (90,270): ax.axvline(x,color="0.75",lw=.8,ls=":")
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11); ax.set_ylabel("signed LOS / sky-mean speed  (+ = receding)",fontsize=11)
ax.set_title("Expanding vs free-fall collapse (red = receding, blue = approaching)",fontsize=12)
ax.legend(frameon=False,fontsize=10,ncol=1,loc="upper center",bbox_to_anchor=(0.5,-0.13)); f.savefig("ff_vs_exp_cc.png",dpi=130,bbox_inches="tight")
print("exp",desc(Se))
