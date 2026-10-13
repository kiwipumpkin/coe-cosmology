"""Collapse level sets. Each dust shell falls from its own stall radius r_max(r) = kappa*D*(r/D)^m under enclosed mass M(r):
 v(r) = -sqrt(2 G M(r) (1/r - 1/r_max(r)))   (exact energy conservation for non-crossing shells; v=0 if shell is at/just past stall: r_max<=r)
kappa = observer's stall radius / D; m = how stall radius scales with shell radius (m=1: all shells same r_max/r, i.e. same shape as r_max=inf).
Mass: point mass 7.5e52 kg (+ optional uniform shell rho=2.9e-27 between 43.9 and 72.5 Gpc, scaled with D).
Fractions solid-angle weighted; Chris convention 0deg = away from CoE."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
G=6.674e-11; GPC=3.0857e25; MP=7.52e52; RHO=2.9e-27
th=np.radians(np.linspace(0,180,721)); ct,st=np.cos(th),np.sin(th); w=st/st.sum()
def Menc(r,D,shell):
    if not shell: return MP+0*r
    a,b=43.9*D/58.2,72.5*D/58.2; rr=np.clip(r,a,b); return MP+4/3*np.pi*RHO*((rr*GPC)**3-(a*GPC)**3)
def vfun(D,kappa,m,shell):
    def v(r):
        r=np.asarray(r,float); rmax=kappa*D*(r/D)**m
        q=2*G*Menc(r,D,shell)*(1/(r*GPC)-1/(rmax*GPC)); return -np.sqrt(np.maximum(q,0))
    return v
def vlos(vf,D,s,c=ct,sn=st):
    x=D+s*c; y=s*sn; r=np.hypot(x,y); return vf(r)*(x*c+y*sn)/r-vf(np.array(D))*c
def stats(D,s,kappa,m,shell=False):
    v=vlos(vfun(D,kappa,m,shell),D,s); return float((w*(v>0)).sum()),float((w*v).sum()/1e3)
# sweep 1: kappa vs m (D=58.2, s=14.3, point mass)
K=np.linspace(1.0,4.0,61); Mm=np.linspace(-4,3,71)
def sweep(xs,ys,fn):
    R=np.zeros((len(ys),len(xs))); V=np.zeros_like(R)
    for i,y in enumerate(ys):
        for j,x in enumerate(xs): R[i,j],V[i,j]=fn(x,y)
    return R,V
R1,V1=sweep(K,Mm,lambda k,m:stats(58.2,14.3,k,m))
R2,V2=sweep(np.linspace(1,30,59),Mm,lambda s,m:stats(58.2,s,1.05,m))
R3,V3=sweep(K,Mm,lambda k,m:stats(58.2,14.3,k,m,shell=True))
for nm,R,V,xs in (("kappa-m pm",R1,V1,K),("s-m k=1.05",R2,V2,np.linspace(1,30,59)),("kappa-m shell",R3,V3,K)):
    i=np.unravel_index(R.argmax(),R.shape); print(nm,"max receding %.3f at x=%.3g m=%.2f, sky-mean v=%.0f km/s"%(R[i],xs[i[1]],Mm[i[0]],V[i]))
def ls(xs,ys,R,V,xl,title,fn):
    f,ax=plt.subplots(figsize=(8.5,5.6))
    lv=[l for l in np.arange(0.1,1.0,0.1) if R.min()<l<R.max()]
    cs=ax.contour(xs,ys,R,levels=lv,colors="#c0392b",linewidths=1.1); ax.clabel(cs,fmt=lambda v:f"{v*100:.0f}%",fontsize=9)
    if R.min()<.5<R.max(): c5=ax.contour(xs,ys,R,levels=[.5],colors="#c0392b",linewidths=3); ax.clabel(c5,fmt="50%%",fontsize=10)
    vm=np.nanmax(np.abs(V)); vl=[x for x in np.linspace(-vm,vm,11) if abs(x)>1e-9]
    cv=ax.contour(xs,ys,V,levels=vl,colors="#2e5fa8",linewidths=1,linestyles="--"); ax.clabel(cv,fmt=lambda v:f"{v:,.0f}",fontsize=8)
    if V.min()<0<V.max(): c0=ax.contour(xs,ys,V,levels=[0],colors="#2e5fa8",linewidths=2.6,linestyles="--"); ax.clabel(c0,fmt="mean 0",fontsize=9)
    ax.set_xlabel(xl,fontsize=11); ax.set_ylabel("stall-radius exponent m  (r_max ∝ r^m)",fontsize=11); ax.set_title(title,fontsize=11.5); ax.grid(alpha=.25)
    ax.plot([],[],color="#c0392b",label="receding (redshift) sky fraction"); ax.plot([],[],color="#2e5fa8",ls="--",label="sky-mean LOS velocity (km/s)")
    ax.legend(frameon=False,fontsize=9.5,loc="upper center",bbox_to_anchor=(0.5,-0.12),ncol=2); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
ls(K,Mm,R1,V1,"observer stall radius / D  (κ)","Collapse, point mass 7.5×10⁵² kg, D = 58.2, s = 14.3 Gpc","collapse_ls_kappa_m.png")
ls(np.linspace(1,30,59),Mm,R2,V2,"galaxy distance s (Gpc)","Collapse, point mass, κ = 1.05 (observer near its stall), D = 58.2","collapse_ls_s_m.png")
ls(K,Mm,R3,V3,"observer stall radius / D  (κ)","Collapse, point mass + uniform shell (2.9×10⁻²⁷ kg/m³), D = 58.2, s = 14.3","collapse_ls_kappa_m_shell.png")
# example curves 0-360
t=np.radians(np.linspace(0,360,1441)); c2,s2=np.cos(t),np.sin(t); w2=np.abs(s2)/np.abs(s2).sum()
ex=[]
for nm,R,xs,sh,ss in (("pm",R1,K,False,None),("shell",R3,K,True,None)):
    i=np.unravel_index(R.argmax(),R.shape); ex.append((xs[i[1]],Mm[i[0]],sh,14.3))
ex.append((4.0,1.0,False,14.3))  # reference: same shape as r_max = inf
def curves(exl,fn,title):
    f,ax=plt.subplots(figsize=(9,5)); ax.axhline(0,color="k",lw=1.2); cols=["#c0392b","#8e44ad","#2e5fa8"]
    for (k,m,sh,s),cl in zip(exl,cols):
        v=vlos(vfun(58.2,k,m,sh),58.2,s,c2,s2)/1e3; fr=(w2*(v>0)).sum(); mean=(w2*v).sum()
        lab=f"κ = {k:.2f}, m = {m:.2f}, {'point mass + shell' if sh else 'point mass'}: receding {fr*100:.0f}% of sky, sky-mean {mean:,.0f} km/s"
        ax.plot(np.degrees(t),v,color=cl,lw=2.2,label=lab); ax.axhline(mean,color=cl,lw=1.2,ls="--")
        print(lab)
    lo,hi=ax.get_ylim(); ax.axhspan(0,hi,color="#fbe3e3",zorder=0); ax.axhspan(lo,0,color="#e3ebfb",zorder=0); ax.set_ylim(lo,hi)
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11)
    ax.set_ylabel("LOS relative velocity (km/s, + = receding)",fontsize=11); ax.set_title(title,fontsize=12)
    ax.legend(frameon=False,fontsize=9,loc="upper center",bbox_to_anchor=(0.5,-0.13)); f.savefig(fn,dpi=130,bbox_inches="tight")
curves([ex[0],(1.05,0.0,False,14.3),ex[2]],"collapse_ls_curves_360.png","Collapse, point mass: smallest blueshift windows (dashed = sky mean)")
curves([ex[1]],"collapse_ls_curves_shell_360.png","Collapse, point mass + uniform shell (dashed = sky mean)")
