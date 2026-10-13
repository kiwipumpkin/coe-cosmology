import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("integ.py").read().split("Ms=[0,1,2,4,8]")[0])
def field(sg,m,kk,TT):
    ri=np.linspace(2,300,3000)*GPC; GM=m*GM1; r=ri.copy(); v=sg*kk*ri; N=4000; h=TT/N
    f=lambda r:-GM/np.maximum(r,1e20)**2
    for _ in range(N):
        k1r,k1v=v,f(r); k2r,k2v=v+.5*h*k1v,f(r+.5*h*k1r); k3r,k3v=v+.5*h*k2v,f(r+.5*h*k2r); k4r,k4v=v+h*k3v,f(r+h*k3r)
        r=r+h/6*(k1r+2*k2r+2*k3r+k4r); v=v+h/6*(k1v+2*k2v+2*k3v+k4v)
    rf=r/GPC; mono=np.all(np.diff(rf)>0) and rf.min()<D-R and rf.max()>D+R
    return (lambda x:np.interp(x,rf,v)), GM, mono
rin=np.linspace(D-R,D+R,400)
def calc(m,kk,TT):
    vf,GM,mono=field(-1,m,kk,TT); vr=vf(rp)[:,None]*xp/rp[:,None]-vf(np.array(D))*xo/D
    sv=vf(rin)/c; return (vr*n).sum(1)/c, np.linalg.norm(vr,axis=1), vf(np.array(D))/c, sv.max(), np.abs(sv).max(), mono
mean=lambda y:(w*y).sum()/w.sum(); rec=lambda y:100*(w*(y>0)).sum()/w.sum()
import sys
if sys.argv[1:]==["scan"]:
    for TT in (174,87,40):
        for fac in (0.25,0.5,1):
            for m in (0,10):
                vl,_,vus,smax,sabs,mono=calc(m,fac*k,TT*GYR); print("T%g k%.2f M%g: v_us %+.3f max shell v %+.3f max|v| %.3f mono %s minmax vLOS %.3f..%.3f"%(TT,fac,m,vus,smax,sabs,mono,vl.min(),vl.max()))
TT=T
def fam(items,title,fv,fa):
    cols=plt.cm.viridis(np.linspace(0.05,0.92,len(items)-1)); data=[]
    for lab,m,kk in items:
        vl,va,vus,smax,sabs,mono=calc(m,kk,TT); v0=calc(0,kk,TT)[1]
        assert mono and smax<0 and sabs<1, (lab,mono,smax,sabs)
        A=va/v0; data.append((lab,m,vl,A,vus))
        print("%s | %s | v_us %+.3f | max shell v %+.3f max|v| %.3f | vLOS %.4f..%.4f mean %+.4f rec %.1f%% | A %.3f..%.3f mean %.4f"%(fv,lab,vus,smax,sabs,vl.min(),vl.max(),mean(vl),rec(vl),A.min(),A.max(),mean(A)))
    for idx,fn,yl,yr in ((2,fv,"v_LOS / c  (+ = receding)",(-0.2,0.33)),(3,fa,"A / A_all = |v_rel| / |v_rel|(M=0, same k)",(0,2.0))):
        f,ax=plt.subplots(figsize=(9,5.3)); ax.axhline(1 if idx==3 else 0,color="0.6",lw=.8); j=0
        for d in data:
            y=d[idx]; L=d[0]+" | v_us %+.3f c | mean %+.4g"%(d[4],mean(y))+(" c, rec %.0f%%"%rec(y) if idx==2 else "")
            if d[1]==0: ax.plot(ang,y,"k--",lw=1.6,label=L)
            else: ax.plot(ang,y,color=cols[j],lw=2,label=L); j+=1
        lo=min(min(d[idx].min() for d in data),yr[0]); hi=max(max(d[idx].max() for d in data),yr[1]); ax.set_ylim(lo-.02*(hi-lo),hi+.02*(hi-lo))
        ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE"); ax.set_ylabel(yl)
        ax.set_title(title,fontsize=10); ax.legend(frameon=False,fontsize=8.5,loc="upper center",bbox_to_anchor=(.5,-.13))
        f.savefig(fn,dpi=125,bbox_inches="tight"); plt.close(f)
if sys.argv[1:]!=["scan"]:
    kc=lambda kk:kk*GPC/c; K=0.5*k
    fam([("no mass",0,K)]+[("M = %g×"%m,m,K) for m in (1,2,4,8,10)],"COLLAPSING — integrated shells, start v = −k r, k = 0.5 k0 = %.4f c/Gpc, T = %.0f Gyr; all shells in view inward, sub-c\nobserver at D = 58.2, R = 14.3 Gpc at final time; M in 7.5×10⁵² kg"%(kc(K),TT/GYR),"integ_vlos_collapse_mass_v2.png","integ_aaall_collapse_mass_v2.png")
    fam([("no mass, k = 0.5 k0",0,K)]+[("k = %g k0 = %.4f c/Gpc"%(f_,kc(f_*k)),10,f_*k) for f_ in (0.25,0.5,0.75)],"COLLAPSING — integrated shells, M = 10× fixed, level sets of start k (k0 = %.4f c/Gpc); all shells in view inward, sub-c\nT = %.0f Gyr, observer at D = 58.2, R = 14.3 Gpc"%(kc(k),TT/GYR),"integ_vlos_collapse_kset_v2.png","integ_aaall_collapse_kset_v2.png")
