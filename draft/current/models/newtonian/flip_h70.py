import sys,numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open('levelD_h70.py').read().split("plt.rcParams")[0].replace("DS=[0.8,1.6,2.4,3.2]","DS=[3.2]"))
KT=float(sys.argv[1]); Ms=[0]+[float(x) for x in sys.argv[2].split(",")]; PLOT=len(sys.argv)>3
MODES["collapse"]=(H*(1-KT),KT/(H*(1-KT)),-1,np.linspace(0.05,20,8000))
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":14,"xtick.labelsize":16,"ytick.labelsize":16})
D=3.2; RM=0.246*D; xo=np.array([D,0.]); Rs=np.linspace(0.01,RM,150)
ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
print("kT=%.2f k=%.1f km/s/Mpc T=%.2f Gyr"%(KT,H*(1-KT)*GPC/c*KMSMPC,KT/(H*(1-KT))/3.156e16))
f2,ax2=plt.subplots(figsize=(8,6)); cc=plt.cm.viridis(np.linspace(.1,.9,len(Ms)-1)); j=0
for m in Ms:
    rf,v=field("collapse",m); sel=(rf>D-RM-.01)&(rf<D+RM+.01); i1=np.where(sel)[0].max()+1
    print("M=%g (%.2e kg): inward %s, no-cross %s, v %.3f..%.3fc, our %.3fc"%(m,m*Mscale*7.5e52,np.all(v[sel]<0),np.all(np.diff(rf[:i1+1])>0),v[sel].min()/c,v[sel].max()/c,np.interp(D,rf,v)/c))
    Fm=lambda x:np.interp(x,rf,v)
    def vl(R,nn):
        xp=xo+np.multiply.outer(R,nn) if nn.ndim==1 else xo+R*nn; rp=np.linalg.norm(xp,axis=-1)
        return ((Fm(rp)[...,None]*xp/rp[...,None]-Fm(D)*xo/D)*nn).sum(-1)/1e3
    y=(vl(RM+.002,n)-vl(RM-.002,n))/4; mn=(w*y).sum()/w.sum(); sd=np.sqrt((w*(y-mn)**2).sum()/w.sum())
    print("   H 0° %.1f 90° %.1f 180° %.1f mean %.1f sd %.1f pos %.1f%%"%(y[0],y[360],y[720],mn,sd,100*(w*(y>0)).sum()/w.sum()))
    if not PLOT: continue
    if m==0: ax2.plot(ang,y,"k--",lw=2.5,label="no mass")
    else: ax2.plot(ang,y,color=cc[j],lw=3.5,label="%g×"%m); j+=1
    f,ax=plt.subplots(figsize=(8,6))
    for a,cl in zip([0,45,90,135,180],["#b2182b","#ef8a62","#4d4d4d","#67a9cf","#2166ac"]):
        ax.plot(Rs*1000,vl(Rs,np.array([np.cos(np.radians(a)),np.sin(np.radians(a))])),color=cl,lw=3.5,label="%d°"%a)
    ax.plot(Rs*1000,-70*Rs*1000,"k--",lw=2,label="no mass"); ax.axhline(0,color="gray",lw=1)
    ax.set_xlabel("galaxy distance R (Mpc)"); ax.set_ylabel("v_LOS (km/s)"); ax.set_xlim(0,RM*1000)
    ax.set_title("Collapsing, "+("no mass" if m==0 else "M = %g×"%m)); ax.legend(frameon=False,ncol=2,fontsize=13)
    f.tight_layout(); f.savefig("h70flip_collapse_angles_v_M%g.png"%m,dpi=110); plt.close(f)
if PLOT:
    ax2.axhline(0,color="gray",lw=1); ax2.set_xlim(0,360); ax2.set_xticks(range(0,361,90)); ax2.set_xlabel("angle (deg), 0° = away from CoE"); ax2.set_ylabel("dv/dR (km/s/Mpc)")
    ax2.set_title("Collapsing, R = %.0f Mpc"%(RM*1000)); ax2.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18),fontsize=13)
    f2.tight_layout(); f2.savefig("h70flip_collapse_slope_vs_angle.png",dpi=110)
