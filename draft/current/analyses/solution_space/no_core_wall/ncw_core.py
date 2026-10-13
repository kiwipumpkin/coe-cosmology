"""No hard CoE-interior wall. Model = ../ss_core.Model (unchanged) with violations() replaced.
Density (the model's own, local_web/web_core.rho_gpc):  rho(r) = RHO * exp(-(r - D)/L),  r = distance from the CoE,
RHO = 2.9e-27 kg/m^3, D = 58.2 Gpc, L = 922.65 Gpc. Smooth everywhere, including inside 43.9 Gpc (rho(0) = 1.065 RHO).
Variants (env NCW = smooth | open):
  open   : the CoE-interior exclusion is simply removed. Everything else as before (cores outside our 14.3 Gpc sphere,
           cores don't overlap, mass 1e49..2e52 kg, void radius 0.5..4 Gpc, centres <= 43 Gpc from us).
  smooth : mass limited by what the local smooth density can supply:
           * cap  M <= 2e52 kg * rho(x)/rho(D)   (the old cap, scaled with the local density)
           * supply volume: a mass M gathers its material from a sphere V_s = M / rho(x) (R_s = (3M/4pi rho)^(1/3),
             ~3.8 Gpc for 2e52 kg). Supply spheres and void spheres may not overlap each other (no material counted twice)
             and may not reach into our visible sphere.  Void deficit = 0.7 rho(x) V already follows the local density.
Net-mass compensation (FC pattern) is kept exactly as in ss_core for comparability."""
import os,sys,numpy as np
ME=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(ME))
from ss_core import *
VAR=os.environ.get("NCW","smooth")
RHO_D=RHO
class NCW(Model):
    def supply(self,p):
        X,M,Rv,Rc=self.decode(p); isC=self.k=="C"
        rho=rho_gpc(X); Rs=np.where(isC,(3*np.abs(M)/(4*np.pi*rho*GPC**3))**(1/3),Rv)
        return X,M,Rv,Rc,rho,Rs,isC
    def violations(self,p):
        if VAR=="open":
            X,M,Rv,Rc=self.decode(p); v=[np.zeros(1)]
            if self.n>1:
                i,j=np.triu_indices(self.n,1); v.append(np.maximum(0,Rc[i]+Rc[j]+GAP-np.linalg.norm(X[i]-X[j],axis=1)))
            return np.concatenate(v)
        X,M,Rv,Rc,rho,Rs,isC=self.supply(p)
        v=[np.where(isC,np.maximum(0,np.log10(M/(MHI*rho/RHO_D))),0.0)*10,   # density-scaled cap (0.1 dex -> 1)
           np.maximum(0,RMIN+Rs-np.linalg.norm(X,axis=1))]                  # supply/void sphere outside our sphere
        if self.n>1:
            i,j=np.triu_indices(self.n,1); v.append(np.maximum(0,Rs[i]+Rs[j]+GAP-np.linalg.norm(X[i]-X[j],axis=1)))
        return np.concatenate(v)
def core_depth(X):
    """how far inside the old wall (43.9 Gpc from CoE) a point is (Gpc, >0 = inside)"""
    return A_IN_G-np.linalg.norm(np.atleast_2d(X)-C_G,axis=1)
