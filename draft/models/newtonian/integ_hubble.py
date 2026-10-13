import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("integ_hik.py").read().split("def calc")[0])     # field(sg,m,kk) RK4, k, T, D, ang, n, w
K=3*k; Ms=[0,1,4,10]; F={m:field(1,m,K)[0] for m in Ms}
xo=np.array([D,0.]); KMSMPC=299792.458/1000
def vlos(m,Rv,nn):
    xp=xo+Rv*nn; rp=np.linalg.norm(xp,axis=-1); vf=F[m]
    vr=vf(rp)[...,None]*xp/rp[...,None]-vf(np.array(D))*xo/D; return (vr*nn).sum(-1)/c
def Rmax(m):
    for Rv in np.arange(14.3,40,0.5):
        if F[m](np.linspace(D-Rv,D+Rv,300)).min()<=0: return Rv-0.5
    return 40
Rlim={m:min(Rmax(m),25) for m in Ms}; print("R limit (outward filter, capped 25):",Rlim)
Rs=np.linspace(0.5,25,200); dirs=[0,45,90,135,180]; cols=plt.cm.coolwarm(np.linspace(1,0,5))
f,axs=plt.subplots(1,4,figsize=(18,4.6),sharey=True)
for ax,m in zip(axs,Ms):
    R_=Rs[Rs<=Rlim[m]]
    for a_,cl in zip(dirs,cols):
        nn=np.array([np.cos(np.radians(a_)),np.sin(np.radians(a_))]); ax.plot(R_,[vlos(m,r,nn) for r in R_],color=cl,lw=2,label="θ = %d°"%a_)
    ax.axvline(14.3,color="0.6",ls=":"); ax.axhline(0,color="k",lw=.7); ax.set_xlabel("galaxy distance R (Gpc)")
    ax.set_title(("no mass" if m==0 else "M = %g× 7.5×10⁵² kg"%m)+"\n(shown while all shells in view outward)",fontsize=9)
axs[0].set_ylabel("v_LOS / c"); axs[0].legend(frameon=False,fontsize=8)
f.suptitle("EXPANDING (integrated, k = 3 k0, T = 174 Gyr, observer at 58.2 Gpc): Hubble diagram by direction (0° = away from CoE)",fontsize=11)
f.tight_layout(); f.savefig("integ_hubble_expand_diagram.png",dpi=115); plt.close(f)
# slopes
nn=n; h=0.05
def slope(m,Rv): return (vlos(m,Rv+h,nn)-vlos(m,Rv-h,nn))/(2*h)
mean=lambda y:(w*y).sum()/w.sum()
f,axs=plt.subplots(1,2,figsize=(16,5),sharey=True); cm=plt.cm.viridis(np.linspace(.1,.9,3))
print("slope table: c/Gpc (km/s/Mpc)")
for ax,Rv in zip(axs,(1.0,14.3)):
    j=0
    for m in Ms:
        s=slope(m,Rv); mn=mean(s); an=100*(s.max()-s.min())/mn
        print("R=%g M=%g: 0° %.5f (%.1f) 90° %.5f (%.1f) 180° %.5f (%.1f) | sky mean %.5f (%.1f) | max %.1f@%.0f min %.1f@%.0f | anisotropy (max-min)/mean %.1f%%"%(Rv,m,s[0],s[0]*KMSMPC,s[360],s[360]*KMSMPC,s[720],s[720]*KMSMPC,mn,mn*KMSMPC,s.max()*KMSMPC,ang[s.argmax()],s.min()*KMSMPC,ang[s[:721].argmin()],an))
        y=s*KMSMPC; L=("no mass" if m==0 else "M = %g×"%m)+" | sky mean %.1f km/s/Mpc | anisotropy %.0f%%"%(mn*KMSMPC,an)
        if m==0: ax.plot(ang,y,"k--",lw=1.6,label=L)
        else: ax.plot(ang,y,color=cm[j],lw=2,label=L); j+=1
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE")
    ax.set_title("local Hubble rate dv_LOS/dR at R = %g Gpc"%Rv,fontsize=10); ax.legend(frameon=False,fontsize=8.5,loc="upper center",bbox_to_anchor=(.5,-.13))
    sec=ax.secondary_yaxis("right",functions=(lambda y:y/KMSMPC,lambda y:y*KMSMPC)); sec.set_ylabel("c per Gpc")
axs[0].set_ylabel("dv_LOS/dR  (km/s/Mpc)")
f.suptitle("EXPANDING (integrated, k = 3 k0, T = 174 Gyr): effective Hubble rate vs angle, mass level sets",fontsize=11)
f.tight_layout(); f.savefig("integ_hubble_expand_slope.png",dpi=115)
