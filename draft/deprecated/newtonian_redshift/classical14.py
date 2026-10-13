import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, pickle
exec(open('integ_h70D.py').read().split('ang=')[0])
kT40=199.0; k=H*(1+kT40); T=kT40/k
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":14,"xtick.labelsize":16,"ytick.labelsize":16})
A=np.arange(0,181,1.0); th=np.radians(A); w=np.sin(th); n=np.stack([np.cos(th),np.sin(th)],1)
DS=[20,30,40,50,60]; Ms=[0,1,4,10]; N=20000; h=T/N; RP=14.3; RES={}
out=open("classical14_results.txt","w")
for D in DS:
    Msc=(k/kp)**2*(D/58.2)**3*((1+kTp)/(1+kT40))**3; ri=np.concatenate([[0.0],np.linspace(1e-6,4.0,8000)])*D*GPC/(1+kT40)
    for m in Ms:
        GM=m*Msc*GM1; r=ri.copy(); v=k*ri; f=lambda r:-GM/np.maximum(np.abs(r),1e18)**2
        for i in range(N):
            a1,q1=v,f(r); a2,q2=v+.5*h*q1,f(r+.5*h*a1); a3,q3=v+.5*h*q2,f(r+.5*h*a2); a4,q4=v+h*q3,f(r+h*a3)
            r=r+h/6*(a1+2*a2+2*a3+a4); v=v+h/6*(q1+2*q2+2*q3+q4)
        mk=(r>0)&(r>=np.maximum.accumulate(r)); R_=r[mk]; V_=v[mk]
        xo=np.array([D*GPC,0.]); uo=np.interp(D*GPC,R_,V_)*np.array([1.,0.])
        xe=xo+RP*GPC*n; re=np.linalg.norm(xe,axis=1); ue=np.interp(re,R_,V_)[:,None]*xe/re[:,None]
        vlos=((ue-uo)*n).sum(1); y=vlos/(RP*1e3)/1e3  # m/s per Gpc -> km/s/Mpc
        RES[(D,m)]=y; mn=(w*y).sum()/w.sum(); an=100*(y.max()-y.min())/mn
        s="D=%d M=%gx  0° %.1f 45° %.1f 90° %.1f 135° %.1f 180° %.1f mean %.1f anis %.0f%%  (emitter v_los max %.2fc)"%(D,m,y[0],y[45],y[90],y[135],y[180],mn,an,vlos.max()/c); print(s); out.write(s+"\n")
out.close(); pickle.dump(RES,open("classical14_res.pkl","wb"))
CZ=pickle.load(open("cosmoz14_res.pkl","rb")); Rg=np.linspace(0.05,15.5,310); ip=np.argmin(abs(Rg-RP))
cc=plt.cm.viridis(np.linspace(.1,.9,3))
def leg(ax): ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18))
def fin(ax,t,f): ax.set_xlim(0,180); ax.set_xticks(range(0,181,45)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_title(t); leg(ax); ax.figure.tight_layout(); ax.figure.savefig(f,dpi=110); plt.close(ax.figure)
for D in DS:
    fig,ax=plt.subplots(figsize=(8,6)); j=0
    for m in Ms:
        y=RES[(D,m)]
        if m==0: ax.plot(A,y,"k--",lw=2.5,label="no mass")
        else: ax.plot(A,y,color=cc[j],lw=3.5,label="%g×"%m); j+=1
    ax.set_ylabel("H = v_los/R (km/s/Mpc)"); fin(ax,"Classical, D = %d Gpc, R = 14.3 Gpc"%D,"classical_angle_D%d.png"%D)
for D in (20,60):
    fig,ax=plt.subplots(figsize=(8,6)); j=0
    ax.plot(A,RES[(D,0)],"k--",lw=2.5,label="no mass")
    for m in Ms[1:]:
        ax.plot(A,CZ[(D,m)][0][ip],color=cc[j],lw=3.5,label="%g× cosmo"%m); ax.plot(A,RES[(D,m)],color=cc[j],lw=3,ls=":",label="%g× classical"%m); j+=1
    ax.set_ylabel("H (km/s/Mpc)"); ax.legend(frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(.5,-.18),fontsize=13)
    ax.set_xlim(0,180); ax.set_xticks(range(0,181,45)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_title("D = %d Gpc: solid cosmological, dotted classical"%D,fontsize=16); fig.tight_layout(); fig.savefig("compare_D%d.png"%D,dpi=110); plt.close(fig)
for tag,yl in (("mean","sky-mean H (km/s/Mpc)"),("anis","anisotropy (max−min)/mean, %")):
    fig,ax=plt.subplots(figsize=(8,6)); j=0
    for m in Ms:
        y=[((w*RES[(D,m)]).sum()/w.sum()) if tag=="mean" else 100*np.ptp(RES[(D,m)])/((w*RES[(D,m)]).sum()/w.sum()) for D in DS]
        if m==0: ax.plot(DS,y,"k--",lw=2.5,label="no mass")
        else: ax.plot(DS,y,"o-",color=cc[j],lw=3.5,ms=9,label="%g×"%m); j+=1
    ax.set_xlabel("observer distance from CoE D (Gpc)"); ax.set_ylabel(yl); ax.set_title("Classical, R = 14.3 Gpc"); leg(ax); fig.tight_layout(); fig.savefig("classical_%s_vs_D.png"%tag,dpi=110); plt.close(fig)
