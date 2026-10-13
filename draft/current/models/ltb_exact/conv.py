from ltbx import *; import sys
D=float(sys.argv[1]); rs=float(sys.argv[2]); rc=float(sys.argv[3]); NQ=int(sys.argv[4]); NT=int(sys.argv[5]); sub=int(sys.argv[6])
A=np.array([0.,90.,170.,180.]); M=Model(D,rs,rc=rc,NQ=NQ,NT=NT); RL=[6.,8.,10.,12.,14.3]
H,te,re,minr=M.trace(A,RL,sub=sub)
print("D=%g rs=%g rc=%g NQ=%d NT=%d sub=%d"%(D,rs,rc,NQ,NT,sub)," ".join("R=%g:%s"%(R,np.round(H[R],1)) for R in RL)," minr180(label)=%.3f"%minr[-1],flush=True)
