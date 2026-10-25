import time,numpy as np;from core import *
t=time.time();o=np.load('../runs_poly/f20000_p_-0.041_-0.035_0.1353_0.npz')['o']
sm,gm,dist=terms(o,shape=shape_poly([-0.041,-0.035,0.1353,0],2.937))
a1,a2=al0(sm),al0(gm);print(time.time()-t,dist.min(),dist.max())
print(np.round(Dl_al0(a1)[2:7],3),np.round(Dl_map(sm)[2:7],3));print(np.round(Dl_al0(a2)[2:7],3),np.round(Dl_map(gm)[2:7],3));print(np.round(Dl_al0(a1+a2)[2:7],3),np.round(Dl_map(sm+gm)[2:7],3))
