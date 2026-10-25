import sys; sys.path.insert(0,'.')
import numpy as np, grcore as G, grcore_tight as T, iso
from grcore import relsum
mu=0.1; v0=0.99; to=106.0
labs=np.radians(np.r_[0,30,60,np.arange(85,100.01,0.5),120,150,170,180])
def schw(mod,labs):
    ems=[]
    for t in labs:
        e=np.array([np.cos(t),np.sin(t)]); at,h=mod.geodesic(mu,e,relsum(v0,0,e),0.0,to); ems.append(at)
    oat,_=mod.geodesic(mu,np.zeros(2),np.array([v0,0]),0.0,to); Xo,Uo=oat(to)
    return mod.observe(mu,Xo,Uo,ems,to,0.0)[0]
a=schw(G,labs); b=schw(T,labs)
# isotropic-coordinate independent build
Mt=iso.Metric(mu)
def toiso(X):  # Schwarzschild-coordinate Cartesian (mass at (-2,0)) -> isotropic Cartesian about mass
    d=X-np.array([-2.0,0]); r=np.hypot(*d); rho=(r-mu+np.sqrt(r*r-2*mu*r))/2; return d/r*rho
ems=[]
for t in labs:
    e=np.array([np.cos(t),np.sin(t)]); at,h=iso.worldline(Mt,toiso(e),relsum(v0,0,e),0.0,to); ems.append((at,0.0))
oat,_=iso.worldline(Mt,toiso(np.zeros(2)),np.array([v0,0]),0.0,to); Xo,vo=oat(to)
c,info=iso.observe(Mt,Xo,vo,ems,to)
print("deg   standard   tight     isotropic   tight/std-1  iso/std-1")
for d,x,y,z in zip(np.degrees(labs),a,b,c): print("%6.1f %.8f %.8f %.8f %9.1e %9.1e"%(d,x,y,z,y/x-1,z/x-1),flush=True)
np.save('verify99.npy',np.c_[np.degrees(labs),a,b,c])
