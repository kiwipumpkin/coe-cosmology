import sys; sys.path.insert(0,'/workspace/clump_scan/ltbx')
from ltbx import *; import pickle
D=float(sys.argv[1]); rs=float(sys.argv[2]); A=np.arange(0,181,1.)
M=Model(D,rs,rc=0.1,NQ=6000,NT=3000); RL=[1.,2.,3.,4.5,6.,8.,10.,12.,14.3]
H=M.trace2(A,RL+[60.],sub=6); maxd=M.maxd.copy()
res=dict(H={R:H[R] for R in RL},maxd=maxd,edge=M.edge.copy(),cl={R:M.classical(A,R) for R in RL},loc=M.local(A),
 dprop=np.trapezoid(np.interp(np.linspace(M.r[0],M.ro,4000),M.r,M.Rp[-1]/np.sqrt(1+M.E2)),np.linspace(M.r[0],M.ro,4000)),
 vobs=np.interp(D,M.R[-1],M.Rd[-1])/C, Mfrac=float(np.interp(M.ro,M.r,M.GM)/max(M.GMc,1e-300)) if rs>0 else 1.0)
res['HvsR']=M.trace2(np.array([0.,90.,180.]),list(np.arange(0.5,14.31,0.5))+[14.3],sub=6)
pickle.dump(res,open("/workspace/clump_scan/ltbx/sweep/s_D%g_rs%g.pkl"%(D,rs),"wb")); print("done",D,rs,flush=True)
