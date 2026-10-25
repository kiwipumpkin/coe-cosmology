import sys,numpy as np,v5 as V
M=float(sys.argv[1]);R=float(sys.argv[2]);a,b=float(sys.argv[3]),float(sys.argv[4])
ps=np.arange(a,b+0.5,2.0);r,wo,t0=V.trace(M,R,ps);np.save(f'h_M{M:.5f}_R{R}_{int(a)}.npy',np.c_[ps,r]);print(wo,t0,np.isfinite(r[:,0]).sum(),flush=True)
