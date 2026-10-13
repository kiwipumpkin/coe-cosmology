import sys; sys.path.insert(0,'/workspace/clump_scan/ltbx')
from ltbx import *; import pickle
D=float(sys.argv[1]); rs=float(sys.argv[2]); A=np.arange(0,181,1.)
M=Model(D,rs,rc=0.1,NQ=6000,NT=3000); RL=[1.,3.,6.,10.,14.3,60.]
H,te,re,minr=M.trace(A,RL,sub=6)
res=dict(H=H,maxd=M.maxd,cl={R:M.classical(A,R) for R in RL[:-1]},loc=M.local(A),
 dprop=np.trapezoid(np.interp(np.linspace(M.r[0],M.ro,4000),M.r,M.Rp[-1]/np.sqrt(1+M.E2)),np.linspace(M.r[0],M.ro,4000)),vobs=np.interp(D,M.R[-1],M.Rd[-1])/C,
 Hr=None)
# H vs R curve at 0/90/180: finer list for a few angles
H2,_,_,_=M.trace(np.array([0.,90.,180.]),list(np.arange(0.5,14.31,0.5))+[14.3],sub=6)
res['HvsR']=H2
pickle.dump(res,open("/workspace/clump_scan/ltbx/sweep/s_D%g_rs%g.pkl"%(D,rs),"wb")); print("done",D,rs,flush=True)
