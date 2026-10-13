import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, pickle
exec(open('integ_h70D.py').read().split('ang=')[0])
kT40=199.0; k=H*(1+kT40); T=kT40/k
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":14,"xtick.labelsize":16,"ytick.labelsize":16})
A=np.arange(0,181,1.0); th=np.radians(A); w=np.sin(th); n=np.stack([np.cos(th),np.sin(th)],1)
DS=[20,30,40,50,60]; Ms=[0,1,4,10]; N=20000; h=T/N; RP=14.3; RES={}; out=open("rel14_results.txt","w")
for D in DS:
    Msc=(k/kp)**2*(D/58.2)**3*((1+kTp)/(1+kT40))**3; ri=np.concatenate([[0.0],np.linspace(1e-6,4.0,8000)])*D*GPC/(1+kT40)
    for m in Ms:
        GM=m*Msc*GM1; r=ri.copy(); v=k*ri; f=lambda r:-GM/np.maximum(np.abs(r),1e18)**2
        for i in range(N):
            a1,q1=v,f(r); a2,q2=v+.5*h*q1,f(r+.5*h*a1); a3,q3=v+.5*h*q2,f(r+.5*h*a2); a4,q4=v+h*q3,f(r+h*a3)
            r=r+h/6*(a1+2*a2+2*a3+a4); v=v+h/6*(q1+2*q2+2*q3+q4)
        mk=(r>0)&(r>=np.maximum.accumulate(r)); R_=r[mk]; V_=v[mk]
        xo=np.array([D*GPC,0.]); eo=np.interp(D*GPC,R_,V_)/c
        xe=xo+RP*GPC*n; re=np.linalg.norm(xe,axis=1); ee=np.interp(re,R_,V_)/c; rh=xe/re[:,None]
        # rapidity eta=v/c; 1+z = u_e.k / u_o.k, photon travels -n
        opz=(np.cosh(ee)+np.sinh(ee)*(rh*n).sum(1))/(np.cosh(eo)+np.sinh(eo)*n[:,0])
        y=c*np.log(opz)/(RP*GPC)*GPC/1e3/1e3/1e3  # km/s/Mpc
        y=c/1e3*np.log(opz)/(RP*1e3)
        grel=np.cosh(ee)*np.cosh(eo)-np.sinh(ee)*np.sinh(eo)*rh[:,0]; brel=np.sqrt(1-1/grel**2)
        RES[(D,m)]=y; mn=(w*y).sum()/w.sum(); an=100*np.ptp(y)/mn
        s="D=%d M=%gx  0° %.1f 45° %.1f 90° %.1f 135° %.1f 180° %.1f mean %.1f anis %.0f%%  rel speed beta 0°/90°/180° %.3f/%.3f/%.3f"%(D,m,y[0],y[45],y[90],y[135],y[180],mn,an,brel[0],brel[90],brel[180]); print(s); out.write(s+"\n")
out.close(); pickle.dump(RES,open("rel14_res.pkl","wb"))
CL=pickle.load(open("classical14_res.pkl","rb")); CZ=pickle.load(open("cosmoz14_res.pkl","rb")); Rg=np.linspace(0.05,15.5,310); ip=np.argmin(abs(Rg-RP))
cc=dict(zip([1,4,10],plt.cm.viridis(np.linspace(.1,.9,3))))
for D in DS:
  for m in (1,4,10):
    fig,ax=plt.subplots(figsize=(8,6)); col=cc[m]
    ax.axhline(70,color="k",ls="--",lw=2.5,label="no mass")
    ax.plot(A,CZ[(D,m)][0][ip],color="tab:blue",lw=3.5,label="cosmological")
    ax.plot(A,CL[(D,m)],color="tab:orange",lw=3.5,ls=":",label="classical")
    ax.plot(A,RES[(D,m)],color="tab:red",lw=3,ls="-.",label="relativistic")
    ax.set_xlim(0,180); ax.set_xticks(range(0,181,45)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_ylabel("H (km/s/Mpc)")
    ax.set_title("D = %d Gpc, %g× mass, R = 14.3 Gpc"%(D,m)); ax.legend(frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(.5,-.18))
    fig.tight_layout(); fig.savefig("overlay_D%d_M%d.png"%(D,m),dpi=110); plt.close(fig)
