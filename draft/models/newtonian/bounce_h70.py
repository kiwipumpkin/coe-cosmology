import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open('levelD_h70.py').read().split("plt.rcParams")[0].replace("DS=[0.8,1.6,2.4,3.2]","DS=[3.2]"))
import sys; Ms=[0]+[float(x) for x in sys.argv[1].split(",")]; PLOT=len(sys.argv)>2
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":14,"xtick.labelsize":16,"ytick.labelsize":16})
D=3.2; RM=0.246*D; xo=np.array([D,0.]); Rs=np.linspace(0.01,RM,150)
ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
f2,ax2=plt.subplots(figsize=(8,6)); cc=plt.cm.viridis(np.linspace(.1,.9,len(Ms)-1)); j=0
for m in Ms:
    rf,v=field("expand",m); v=-v                      # time reverse
    sel=(rf>D-RM-.01)&(rf<D+RM+.01); GM=m*Mscale*GM1
    E=0.5*v**2-GM/(rf*GPC); Eu=E[sel]>0; rs=rf[sel]
    mono=np.all(np.diff(rf[:np.where(sel)[0].max()+2])>0)
    if m>0: print("M=%g (%.2e kg): inward %s, no-cross %s, |v| %.3f..%.3fc; unbound shells: %s; our E/(GM/r)=%.2f"%(m,m*Mscale*7.5e52,np.all(v[sel]<0),mono,-v[sel].max()/c,-v[sel].min()/c,
        "none" if not Eu.any() else ("ALL" if Eu.all() else "r>%.2f Gpc"%rs[Eu].min()),np.interp(D,rf,E)/(GM/(D*GPC))))
    Fm=(lambda rf,v:lambda x:np.interp(x,rf,v))(rf,v)
    if m==0: Fm=lambda x:-k70*x
    k70=H*GPC                                         # v=-70 r reference (m/s per Gpc)
    if m==0: Fm=lambda x:-k70*x
    def vl(R,nn):
        xp=xo+(np.multiply.outer(R,nn) if nn.ndim==1 else R*nn); rp=np.linalg.norm(xp,axis=-1)
        return ((Fm(rp)[...,None]*xp/rp[...,None]-Fm(D)*xo/D)*nn).sum(-1)/1e3
    for Rv in (0.055,RM):
        y=(vl(Rv+.002,n)-vl(Rv-.002,n))/4; mn=(w*y).sum()/w.sum(); sd=np.sqrt((w*(y-mn)**2).sum()/w.sum())
        print("   R=%3.0f: 0° %.1f 90° %.1f 180° %.1f mean %.1f sd %.1f pos %.1f%%"%(Rv*1000,y[0],y[360],y[720],mn,sd,100*(w*(y>0)).sum()/w.sum()))
    if not PLOT: continue
    if m==0: ax2.plot(ang,y,"k--",lw=2.5,label="no mass")
    else: ax2.plot(ang,y,color=cc[j],lw=3.5,label="%g×"%m); j+=1
    if m==0: continue
    f,ax=plt.subplots(figsize=(8,6))
    for a,cl in zip([0,45,90,135,180],["#b2182b","#ef8a62","#4d4d4d","#67a9cf","#2166ac"]):
        ax.plot(Rs*1000,vl(Rs,np.array([np.cos(np.radians(a)),np.sin(np.radians(a))])),color=cl,lw=3.5,label="%d°"%a)
    ax.plot(Rs*1000,-70*Rs*1000,"k--",lw=2,label="no mass"); ax.axhline(0,color="gray",lw=1)
    ax.set_xlabel("galaxy distance R (Mpc)"); ax.set_ylabel("v_LOS (km/s)"); ax.set_xlim(0,RM*1000)
    ax.set_title("Collapsing after turnaround, M = %g×"%m); ax.legend(frameon=False,ncol=2,fontsize=13)
    f.tight_layout(); f.savefig("bounce_v_M%g.png"%m,dpi=110); plt.close(f)
if PLOT:
    ax2.axhline(0,color="gray",lw=1); ax2.set_xlim(0,360); ax2.set_xticks(range(0,361,90)); ax2.set_xlabel("angle (deg), 0° = away from CoE"); ax2.set_ylabel("dv/dR (km/s/Mpc)")
    ax2.set_title("After turnaround, R = %.0f Mpc"%(RM*1000)); ax2.legend(frameon=False,ncol=5,loc="upper center",bbox_to_anchor=(.5,-.18),fontsize=13)
    f2.tight_layout(); f2.savefig("bounce_slope_Rmax.png",dpi=110)
