import sys;sys.path.insert(0,'.');from grcore import *
mu=0.1
for v0 in (0.9,0.99):
  for to in ([10,10.2,10.5,11,12] if v0==0.9 else [100,101,102,105,110]):
    oat,h=geodesic(mu,np.zeros(2),np.array([v0,0]),0,to); Xo,Uo=oat(to)
    so,L=ray(mu,Xo,0.0,to,0.0); print(v0,to,'observer x %.3f'%Xo[0],'ray at t=0 x %.4f (covered if < -1)'%cart(*so.sol(0.0)[:2])[0],flush=True)
