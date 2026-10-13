"""Cosmological redshift along null rays through the time-dependent radial shell flow (H0=70 setup of integ_h70D.py)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open('integ_h70D.py').read().split('ang=')[0])
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":14,"xtick.labelsize":16,"ytick.labelsize":16})
A=np.arange(0,181,1.0); th=np.radians(A); w=np.sin(th); n=np.stack([np.cos(th),np.sin(th)],1)
Rg=np.linspace(0.05,14.3,286)  # Gpc, present-day distance of emitter
DS=[20,30,40,50,60]; Ms=[0,1,4,10]; RES={}; N=3000; h=T/N; S=10
def cone(y,c0):  # 75deg cone mean, axisym weight sin
    m=np.abs(A-c0)<=75; return (w[m]*y[m]).sum()/max(w[m].sum(),1e-9) if c0 not in (0,180) else None
def capmean(y,c0):
    # proper cone around axis direction c0 in 3D (axisymmetric field): sample sphere
    return CM[c0]@y
# cone weights over 3D sphere about direction at angle c0 (in plane), mapping to polar angle from axis
mu=np.random.default_rng(1).normal(size=(200000,3)); mu/=np.linalg.norm(mu,axis=1)[:,None]
CM={}
for c0 in (0,45,90,135,180):
    d=np.array([np.cos(np.radians(c0)),np.sin(np.radians(c0)),0]); s=mu[mu@d>=np.cos(np.radians(75))]
    pa=np.degrees(np.arccos(np.clip(s[:,0],-1,1))); hst=np.bincount(np.clip(np.rint(pa).astype(int),0,180),minlength=181).astype(float); CM[c0]=hst/hst.sum()
out=open("cosmoz_results.txt","w")
def P(*a):
    s=" ".join(str(x) for x in a); print(s); out.write(s+"\n")
for D in DS:
    Msc=(k/kp)**2*(D/58.2)**3; ri=np.linspace(0.0005,1.6,6000)*D*GPC
    for m in Ms:
        GM=m*Msc*GM1; r=ri.copy(); v=k*ri; f=lambda r:-GM/np.maximum(r,1e18)**2
        snR=[r.copy()]; snV=[v.copy()]
        for i in range(N):
            a1,q1=v,f(r); a2,q2=v+.5*h*q1,f(r+.5*h*a1); a3,q3=v+.5*h*q2,f(r+.5*h*a2); a4,q4=v+h*q3,f(r+h*a3)
            r=r+h/6*(a1+2*a2+2*a3+a4); v=v+h/6*(q1+2*q2+2*q3+q4)
            if (i+1)%S==0: snR.append(r.copy()); snV.append(v.copy())
        snR=np.array(snR); snV=np.array(snV); rnow=snR[-1]
        cross=(np.diff(snR,axis=1)<=0).any(); turn=(snV[-1]<0).any()
        # max local relative speed: neighbours within 100 Mpc proper, any epoch
        dR=np.diff(snR,axis=1); gr=np.abs(np.diff(snV,axis=1)/np.where(dR>0,dR,np.nan)); locmax=np.nanmax(gr*0.1*GPC)/c
        xr=snR[:,1:][dR<=0]; crossmax=(xr.max()/GPC if xr.size else 0); tr=rnow[snV[-1]<0]; turnmax=(tr.max()/GPC if tr.size else 0)
        rO=np.interp(D*GPC,rnow,ri)
        xo=np.array([np.interp(rO,ri,rnow),0.])
        # trace rays backward in time
        x=np.tile(xo,(len(A),1)); lnz=np.zeros(len(A)); Dn=[np.zeros(len(A))]; LZ=[lnz.copy()]
        nt=len(snR)-1
        for j in range(nt,0,-1):
            for sub in range(S):
                tf=(j-(sub+.5)/S)  # fractional snapshot index
                a=j-1; fr=tf-a; R_=snR[a]*(1-fr)+snR[a+1]*fr; V_=snV[a]*(1-fr)+snV[a+1]*fr
                r0=np.linalg.norm(x,axis=1); u0=np.interp(r0,R_,V_)[:,None]*x/r0[:,None]
                xm=x+(n*c-u0)*h/2; rr=np.linalg.norm(xm,axis=1); vr=np.interp(rr,R_,V_); dv=np.interp(rr,R_,np.gradient(V_,R_))
                cos=(xm*n).sum(1)/rr; grad=dv*cos**2+vr/rr*(1-cos**2)
                lnz+=grad*h; um=vr[:,None]*xm/rr[:,None]; x=x+(n*c-um)*h
            # emitter's present position: label at radius |x| at snapshot j-1
            rr=np.linalg.norm(x,axis=1); lab=np.interp(rr,snR[j-1],ri); xn=x/rr[:,None]*np.interp(lab,ri,rnow)[:,None]
            Dn.append(np.linalg.norm(xn-xo,axis=1)); LZ.append(lnz.copy())
        Dn=np.array(Dn)/GPC; LZ=np.array(LZ)
        Rmax=Dn[-1].min()
        Hln=np.full((len(Rg),len(A)),np.nan); Hcz=Hln.copy()
        for ia in range(len(A)):
            z=np.interp(Rg,Dn[:,ia],LZ[:,ia],right=np.nan)
            Hln[:,ia]=c/1e6*z/Rg  # km/s/Mpc: c[m/s]*.../Gpc -> (m/s/Gpc)/1e6 /1e3? fix below
            Hcz[:,ia]=c/1e6*np.expm1(z)/Rg
        RES[(D,m)]=(Hln,Hcz,Rmax)
        P("D=%d M=%gx  horizon (max present distance reachable, min over angles)=%.2f Gpc  shell-crossing=%s (out to r=%.2f Gpc)  infall now out to r=%.2f Gpc  max rel speed within 100 Mpc=%.3fc  our speed now=%.2fc"%(D,m,Rmax,cross,crossmax,turnmax,locmax,np.interp(rO,ri,snV[-1])/c))
        for RR in (1.0,3.0,Rmax*0.98 if Rmax<14.3 else 14.3):
            i=np.argmin(abs(Rg-RR)); y=Hln[i]; yc=Hcz[i]
            if np.isnan(y).any(): continue
            mn=(w*y).sum()/w.sum(); an=100*(y.max()-y.min())/mn; mc=(w*yc).sum()/w.sum()
            P("   R=%.2f Gpc  ln-based H: 0° %.1f 45° %.1f 90° %.1f 135° %.1f 180° %.1f mean %.1f anis %.0f%% | cones75 0:%.1f 45:%.1f 90:%.1f 135:%.1f 180:%.1f | cz-based mean %.1f (0° %.1f 180° %.1f)"%(Rg[i],y[0],y[45],y[90],y[135],y[180],mn,an,*(CM[c0]@y for c0 in (0,45,90,135,180)),mc,yc[0],yc[180]))
out.close()
import pickle; pickle.dump({k_:(v_[0],v_[1],v_[2]) for k_,v_ in RES.items()},open("cosmoz_res.pkl","wb"))

cc=plt.cm.viridis(np.linspace(.1,.9,3)); RP=3.0; ip=np.argmin(abs(Rg-RP))
def leg(ax): ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18))
for D in DS:
    fig,ax=plt.subplots(figsize=(8,6)); j=0
    for m in Ms:
        y=RES[(D,m)][0][ip]
        if m==0: ax.plot(A,y,"k--",lw=2.5,label="no mass")
        else: ax.plot(A,y,color=cc[j],lw=3.5,label="%g×"%m); j+=1
    ax.set_xlim(0,180); ax.set_xticks(range(0,181,45)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_ylabel("H = c·ln(1+z)/R (km/s/Mpc)")
    ax.set_title("Cosmological z, D = %d Gpc, R = 3 Gpc"%D); leg(ax); fig.tight_layout(); fig.savefig("cosmoz_angle_D%d.png"%D,dpi=110); plt.close(fig)
    for m in Ms[1:]:
        Hl,_,Rm=RES[(D,m)]; fig,ax=plt.subplots(figsize=(8,6)); col=plt.cm.plasma(np.linspace(0,.85,5))
        for q,a0 in enumerate((0,45,90,135,180)): ax.plot(Rg,Hl[:,a0],color=col[q],lw=3,label="%d°"%a0)
        mn=(Hl*w).sum(1)/w.sum(); ax.plot(Rg,mn,color="grey",lw=5,alpha=.6,label="sky mean")
        ax.plot(Rg,(RES[(D,0)][0]*w).sum(1)/w.sum(),"k--",lw=2.5,label="no mass")
        ax.set_xlim(0,Rm); ax.set_xlabel("present distance R (Gpc)"); ax.set_ylabel("H = c·ln(1+z)/R (km/s/Mpc)")
        ax.set_title("D = %d Gpc, %g× mass"%(D,m)); ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18),fontsize=13)
        fig.tight_layout(); fig.savefig("cosmoz_HvsR_D%d_M%g.png"%(D,m),dpi=110); plt.close(fig)
for tag,yl in (("mean","sky-mean H (km/s/Mpc)"),("anis","anisotropy (max−min)/mean, %")):
    fig,ax=plt.subplots(figsize=(8,6)); j=0
    for m in Ms:
        y=[]
        for D in DS:
            v_=RES[(D,m)][0][ip]; mn=(w*v_).sum()/w.sum(); y.append(mn if tag=="mean" else 100*(v_.max()-v_.min())/mn)
        if m==0: ax.plot(DS,y,"k--",lw=2.5,label="no mass")
        else: ax.plot(DS,y,"o-",color=cc[j],lw=3.5,ms=9,label="%g×"%m); j+=1
    ax.set_xlabel("observer distance from CoE D (Gpc)"); ax.set_ylabel(yl); ax.set_title("Cosmological z, R = 3 Gpc"); leg(ax)
    fig.tight_layout(); fig.savefig("cosmoz_%s_vs_D.png"%tag,dpi=110); plt.close(fig)
