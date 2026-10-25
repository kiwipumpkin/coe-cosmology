# validate FFT pipeline vs direct point-mass sky() sum on a test set of bodies far from the R=2 sphere
import numpy as np;from build import *
G=setup(dict(name='t',L=20.0,N=384,rin=0,rout=5.0,seed=1));rng=np.random.default_rng(5)
P_=rng.uniform(-4.5,4.5,(400,3));r=np.linalg.norm(P_,axis=1);P_=P_[(abs(r-2)>0.8)&(r<5)&(r>0.3)];m=rng.normal(0,1e19,len(P_))
rho=cic(G,P_,m)/G['dx']**3;s=skyfrom(G,accel(G,rho))
Y=P_/RV;out=np.zeros((len(RS),len(V)))
for i,R in enumerate(RS):
    d=Y[:,None,:]-(R/RV)*V[None];ar=np.einsum('bpk,pk->bp',d,V)/np.linalg.norm(d,axis=2)**3-(Y@V.T)/np.linalg.norm(Y,axis=1)[:,None]**3
    o=(m*MSUN/KG)@ar;out[i]=o-o.mean()
for i,R in enumerate(RS[:4]):print(R,'corr',np.corrcoef(s[i],out[i])[0,1],'ratio rms',s[i].std()/out[i].std())
