"""Sky fraction (solid-angle weighted) with LOS relative velocity >0 / <0. Chris convention 0deg = away from CoE.
Expanding: paper run-2 law with mass scaled: v^2 = b1^2 r0^2 + 2 f mu (1/x - 1/r0), x=r/d1, r0=x/3, d1=19.4 Gpc fixed (physical); invalid if v^2<0 anywhere sampled.
Collapsing: free fall v = -sqrt(2GM(1/r - 1/rmax)), rmax = inf, 1.5x or 2x outer shell (72.5 Gpc scaled with D: r_out = D*72.5/58.2)."""
import numpy as np, matplotlib, json; matplotlib.use("Agg"); import matplotlib.pyplot as plt
b1,mu,d1=0.5866,0.09332,19.4
th=np.radians(np.linspace(0,180,721)); ct,st=np.cos(th),np.sin(th); w=st/st.sum()
def vlos(vf,D,s):
    x=D+s*ct; y=s*st; r=np.hypot(x,y); v=vf(r); vo=vf(np.array(D))
    return v*(x*ct+y*st)/r - vo*ct          # (v_g r_hat - v_o x_hat).n, n=(ct,st), +x away from CoE
def vexp(f): 
    def g(r):
        r=np.asarray(r,float); x=r/d1; q=b1**2*(x/3)**2+2*f*mu*(1/x-3/x); return np.where(q>=0,np.sqrt(np.abs(q)),np.nan)
    return g
def vff(rmax): return lambda r:-np.sqrt(np.maximum(1/np.asarray(r,float)-1/rmax,0))
def frac(v): 
    if np.isnan(v).any(): return np.nan,np.nan
    return float((w*(v>0)).sum()),float((w*(v<0)).sum())
S=np.linspace(1,30,59); DD=np.linspace(20,80,61); F=np.logspace(-1,2,61)
out={}
# collapse: fraction vs s,D for each rmax (mass-independent)
for lab,mult in (("inf",None),("1.5",1.5),("2",2.0)):
    Z=np.array([[frac(vlos(vff(np.inf if mult is None else mult*D*72.5/58.2),D,s))[0] for s in S] for D in DD]); out["col_"+lab]=Z
# expanding: fraction vs mass factor and s at D=58.2; and vs s,D at f=1
out["exp_fs"]=np.array([[frac(vlos(vexp(f),58.2,s))[0] for s in S] for f in F])
out["exp_sD"]=np.array([[frac(vlos(vexp(1),D,s))[0] for s in S] for D in DD])
# mass check for free fall: scale only
v1=vlos(lambda r:-np.sqrt(2*1/r),58.2,14.3); v8=vlos(lambda r:-np.sqrt(2*8/r),58.2,14.3); print("ff mass-independent sign:",np.array_equal(np.sign(v1),np.sign(v8)))
for k,Z in out.items():
    i=np.unravel_index(np.nanargmax(Z),Z.shape); j=np.unravel_index(np.nanargmin(Z),Z.shape)
    ax0=F if k=="exp_fs" else DD
    print(k,"max rec %.3f at (%s=%.3g, s=%.1f)"%(Z[i],"f" if k=="exp_fs" else "D",ax0[i[0]],S[i[1]]),"| min rec %.3f at (%.3g, s=%.1f)"%(Z[j],ax0[j[0]],S[j[1]]),"| nan frac %.2f"%np.isnan(Z).mean())
np.savez("sweep_vel.npz",S=S,DD=DD,F=F,**out)
cm=plt.cm.RdBu_r
def heat(Z,y,ylab,title,fn,logy=False):
    f,ax=plt.subplots(figsize=(8,5.2))
    lv=[l for l in np.arange(0.05,1.0,0.05) if np.nanmin(Z)<=l<=np.nanmax(Z)] or [round(float(np.nanmean(Z)),3)]
    cs=ax.contour(S,y,Z,levels=lv,colors="#2e5fa8",linewidths=1.1); ax.clabel(cs,fmt=lambda v:f"{v*100:.0f}%",fontsize=9)
    if np.nanmin(Z)<0.5<np.nanmax(Z): c5=ax.contour(S,y,Z,levels=[0.5],colors="k",linewidths=2.6); ax.clabel(c5,fmt="50%",fontsize=10)
    if np.isnan(Z).any(): ax.contour(S,y,np.isnan(Z).astype(float),levels=[0.5],colors="0.5",linewidths=1.2,linestyles="--")
    if logy: ax.set_yscale("log")
    ax.set_xlabel("galaxy distance s (Gpc)",fontsize=11); ax.set_ylabel(ylab,fontsize=11); ax.set_title(title,fontsize=12)
    ax.grid(alpha=.25); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
heat(out["col_inf"],DD,"observer distance D from CoE (Gpc)","Free-fall collapse (r_max = ∞): receding-sky-fraction contours\n(any mass — mass only scales speeds)","sweep_vel_col_inf.png")
heat(out["col_2"],DD,"observer distance D from CoE (Gpc)","Free-fall collapse (r_max = 2 × outer shell): receding sky fraction\n(any mass)","sweep_vel_col_rmax2.png")
heat(out["exp_fs"],F,"point mass / 7.5×10⁵² kg","Expanding (run-2 energy law, mass scaled), D = 58.2 Gpc\nreceding sky fraction (grey dashed = edge of valid region; beyond it shells are bound)","sweep_vel_exp_mass.png",logy=True)
# example curves
ang=np.degrees(th)
def curve(items,title,fn):
    f,ax=plt.subplots(figsize=(8.5,4.6)); ax.axhline(0,color="k",lw=1)
    for vf,D,s,col,lab in items:
        v=vlos(vf,D,s); v=v/np.abs(v).max(); ax.plot(ang,v,color=col,lw=2.2,label=lab+f"  (receding {frac(vlos(vf,D,s))[0]*100:.0f}% of sky)")
    lo,hi=ax.get_ylim(); ax.axhspan(0,hi,color="#fbe3e3",zorder=0); ax.axhspan(lo,0,color="#e3ebfb",zorder=0); ax.set_ylim(lo,hi)
    ax.set_xlim(0,180); ax.set_xticks(range(0,181,30)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11)
    ax.set_ylabel("LOS relative velocity / max |v|  (+ = receding)",fontsize=11); ax.set_title(title,fontsize=12)
    ax.legend(frameon=False,fontsize=9.5,loc="upper center",bbox_to_anchor=(0.5,-0.14)); f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
json.dump({k:None for k in out},open("/dev/null","w"))
# r_max sweep at D=58.2 (r_max must exceed every galaxy radius sampled -> r_max >= D+s)
RM=np.linspace(1.0,4.0,61)
Z=np.full((len(RM),len(S)),np.nan)
for i,m in enumerate(RM):
    for j,s in enumerate(S):
        rmax=m*72.5
        if rmax>=58.2+s: Z[i,j]=frac(vlos(vff(rmax),58.2,s))[0]
i=np.unravel_index(np.nanargmax(Z),Z.shape); j=np.unravel_index(np.nanargmin(Z),Z.shape)
print("rmax sweep D=58.2: max rec %.3f at rmax=%.2fx s=%.1f; min %.3f at %.2fx s=%.1f"%(Z[i],RM[i[0]],S[i[1]],Z[j],RM[j[0]],S[j[1]]))
for m in (1.0,1.1,1.25,1.5,2,4): 
    for s in (1,7,14.3): 
        if m*72.5>=58.2+s: print(" rmax=%.2fx s=%.1f rec=%.3f"%(m,s,frac(vlos(vff(m*72.5),58.2,s))[0]))
heat(Z,RM,"stall radius r_max / 72.5 Gpc","Free-fall collapse, D = 58.2 Gpc: receding sky fraction vs stall radius\n(grey dashed = galaxy beyond r_max; any mass)","sweep_vel_col_rmax.png")
curve([(vff(1.0*72.5),58.2,14.3,"#c0392b","free fall, r_max = 72.5 Gpc (stall at outer shell), s = 14.3"),
       (vff(np.inf),58.2,14.3,"#2e5fa8","free fall, r_max = ∞, s = 14.3"),
       (vexp(1),58.2,14.3,"#7f3b3b","expanding, run-2 law, 1× mass, s = 14.3")],"Velocity vs angle at the extremes, D = 58.2 Gpc","sweep_vel_curves.png")
