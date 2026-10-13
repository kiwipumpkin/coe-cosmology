from ltbx import *; import sys, pickle
D=float(sys.argv[1]); rs=float(sys.argv[2]); A=np.arange(0,181,1.)
M=Model(D,rs,rc=0.1,NQ=6000,NT=3000); RL=[1.,3.,6.,8.,10.,12.,14.3]
H,te,re,minr=M.trace(A,RL,sub=6)
res=dict(H=H,te=te,re=re,minr=minr,cl={R:M.classical(A,R) for R in RL},loc=M.local(A),
  dprop=np.trapezoid(np.interp(np.linspace(M.r[0],M.ro,4000),M.r,M.Rp[-1]/np.sqrt(1+M.E2)),np.linspace(M.r[0],M.ro,4000)),vobs=np.interp(D,M.R[-1],M.Rd[-1])/C)
pickle.dump(res,open("ltbx_D%g_rs%g.pkl"%(D,rs),"wb"))
w=np.sin(np.radians(A))
print("D=%g rs=%g  proper dist us->CoE now %.2f Gpc, our speed %.2fc"%(D,rs,res['dprop'],res['vobs']))
for R in RL:
    for nm,y in (("cosmo",H[R]),("class",res['cl'][R])):
        ok=~np.isnan(y); mn=(w[ok]*y[ok]).sum()/w[ok].sum()
        f=lambda a:"%.1f(%.3f)"%(y[a],y[a]/mn) if ok[a] else "beyond-CoE"
        print("  R=%4.1f %s H(A/Aall) 0:%s 45:%s 90:%s 135:%s 180:%s mean %.1f  valid 0-%d"%(R,nm,f(0),f(45),f(90),f(135),f(180),mn,A[ok].max()))
    if R==14.3: print("  z at 180/0: %s / %.1f, t_emit(0) %.2f Gyr"%(np.expm1(H[R][180]*R/KM/C) if ok[180] else '-',np.expm1(H[R][0]*R/KM/C),te[R][0]))
