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
# check against integ.py at k0, 20x
print("check k0 20x:",calc(1,20,k)[0][[0,360,720]])
for fac in (1,1.25,1.5,1.75,2,2.5,3):
    vl,_,vus,vf=calc(1,20,fac*k); print("k=%.2f k0 (%.5f c/Gpc): v_us %+.3f min %.4f@%.0f  v90 %.4f  max start speed in field %.2fc"%(fac,fac*k*GPC/c,vus,vl.min(),ang[vl[:721].argmin()],vl[360],fac*k*GPC/c*(D+R)/(1+0)))
