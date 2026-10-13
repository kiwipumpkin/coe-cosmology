import sys; sys.path.insert(0,'/workspace/clump_scan/ltbx'); from ltbx import *; import pickle
D=float(sys.argv[1]); m=float(sys.argv[2]); M=Model(D,m,rc=0.1,NQ=6000,NT=3000)
A=np.array([0.,45.,89.,89.99,90.,90.01,91.,135.])
H=M.trace2(A,[3.,14.3]); S=pickle.load(open("sweep/s_D%g_rs%g.pkl"%(D,m),"rb"))
pickle.dump(H,open("sweep/t90_D%g_rs%g.pkl"%(D,m),"wb"))
print("D=%g m=%g"%(D,m)," | ".join("R=%g new %s old(0,45,89,90,91,135) %s"%(R,np.round(H[R],3),np.round(S['H'][R][[0,45,89,90,91,135]],3)) for R in (3.,14.3)),flush=True)
