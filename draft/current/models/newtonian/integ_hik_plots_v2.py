import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src=open("integ.py").read().split("Ms=[0,1,2,4,8]")[0]
exec(src)
def field(sg,m,kk):   # fast vectorized RK4 replacement (same ODE)
    ri=np.linspace(2,200,3000)*GPC; GM=m*GM1; r=ri.copy(); v=sg*kk*ri; N=4000; h=T/N
    f=lambda r:-GM/np.maximum(r,1e20)**2
    for _ in range(N):
        k1r,k1v=v,f(r); k2r,k2v=v+.5*h*k1v,f(r+.5*h*k1r); k3r,k3v=v+.5*h*k2v,f(r+.5*h*k2r); k4r,k4v=v+h*k3v,f(r+h*k3r)
        r=r+h/6*(k1r+2*k2r+2*k3r+k4r); v=v+h/6*(k1v+2*k2v+2*k3v+k4v)
    rf=r/GPC; assert np.all(np.diff(rf)>0) and rf.min()<D-R and rf.max()>D+R
    return lambda x:np.interp(x,rf,v), GM
def calc(sg,m,kk):
    vf,GM=field(sg,m,kk); vr=vf(rp)[:,None]*xp/rp[:,None]-vf(np.array(D))*xo/D
    return (vr*n).sum(1)/c, np.linalg.norm(vr,axis=1), vf(np.array(D))/c, vf
mean=lambda y:(w*y).sum()/w.sum(); rec=lambda y:100*(w*(y>0)).sum()/w.sum()

def fam(items,title,fn_v,fn_a):
    cols=plt.cm.viridis(np.linspace(0.05,0.92,len(items)-1))
    data=[]
    for lab,m,kk in items:
        vl,va,vus,vf_=calc(1,m,kk); v0=calc(1,0,kk)[1]; ms=vf_(np.linspace(D-R,D+R,400)).min()/c; assert ms>0 and vl.min()>0, (lab,ms,vl.min()); print('  filter OK: min shell v %.4f c, min vLOS %.4f c'%(ms,vl.min())); data.append((lab,m,vl,va/v0,vus))
        print(fn_v,lab,"v_us %+.3f min %.4f@%.0f v90 %.4f mean %+.4f rec %.0f%% | A min %.3f max %.3f mean %.4f"%(vus,vl.min(),ang[vl[:721].argmin()],vl[360],mean(vl),rec(vl),(va/v0).min(),(va/v0).max(),mean(va/v0)))
    for idx,fn,yl in ((2,fn_v,"v_LOS / c  (+ = receding)"),(3,fn_a,"A / A_all = |v_rel| / |v_rel|(M=0, same k)")):
        f,ax=plt.subplots(figsize=(9,5.3)); ax.axhline(1 if idx==3 else 0,color="0.6",lw=.8); j=0
        for d in data:
            y=d[idx]; L=d[0]+" | v_us %+.3f c | mean %+.4g"%(d[4],mean(y))+(" c, rec %.0f%%"%rec(y) if idx==2 else "")
            if d[1]==0: ax.plot(ang,y,"k--",lw=1.6,label=L)
            else: ax.plot(ang,y,color=cols[j],lw=2,label=L); j+=1
        ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE"); ax.set_ylabel(yl)
        ax.set_title(title,fontsize=10); ax.legend(frameon=False,fontsize=8.5,loc="upper center",bbox_to_anchor=(.5,-.13))
        f.savefig(fn,dpi=125,bbox_inches="tight"); plt.close(f)

K3=3*k; kc=lambda kk:kk*GPC/c
fam([("no mass",0,K3)]+[("M = %g×"%m,m,K3) for m in (1,2,4,8,10)],"EXPANDING — integrated shells, start k = 3 k0 = %.4f c/Gpc, T = %.0f Gyr; all shells in view moving outward\nobserver at D = 58.2, R = 14.3 Gpc at final time; M in 7.5×10⁵² kg"%(kc(K3),T/GYR),"integ_vlos_expand_hik_c_v2.png","integ_aaall_expand_hik_c_v2.png")
fam([("no mass, k = 3 k0",0,K3)]+[("k = %g k0 = %.4f c/Gpc"%(f_,kc(f_*k)),10,f_*k) for f_ in (3,4,5,6)],"EXPANDING — integrated shells, M = 10× fixed, level sets of start k (k0 = %.4f c/Gpc); all shells in view outward\nT = %.0f Gyr, observer at D = 58.2, R = 14.3 Gpc"%(kc(k),T/GYR),"integ_vlos_expand_kset_c_v2.png","integ_aaall_expand_kset_c_v2.png")
