"""ALL OPEN: objects anywhere 0.5 Gpc .. 43 Gpc from us -- inside our 14.3 Gpc sphere, in the shell, or in the old CoE core.
Kernel: ss_core.unit_field (= lattice_core.field3) for every object. It is plain Newtonian (pull at ring point - pull on us,
  projected on the sightline, mean removed) and was verified for interior points in inside_visible (corr >= 0.99999999),
  so the 'inside' and 'outside' kernels are the same function. Objects may NOT straddle the ring: each core / void sphere must
  lie wholly inside (d + R <= 14.2) or wholly outside (d - R >= 14.4) -- straddling is forbidden, not split.
Us: no core / void sphere within 0.5 Gpc of us (d >= R + 0.5).
Density: rho(r) = 2.9e-27 kg/m^3 * exp(-(r - 58.2)/922.65), r = distance from CoE (web_core.rho_gpc). No wall anywhere.
  smooth : mass cap M <= 2e52 kg * rho(x)/rho_ref; each mass draws its material from a supply sphere V = M/rho(x);
           supply spheres and void spheres may not overlap each other (no material counted twice).
           (supply spheres may cross the ring and enclose us -- only the physical cores/voids are restricted there.)
  open   : no density limit; usual caps (1e49..2e52 kg, voids 0.5..4 Gpc), cores don't overlap.
Compensation: inside_visible 'local' convention for every fit -- net mass removed smoothly over our visible sphere with the
  exp-lock density shape (inside_core.FLOC, analytic gradient dipole; the uniform part gives exactly zero)."""
import os,sys,numpy as np
ME=os.path.dirname(os.path.abspath(__file__)); SS=os.path.dirname(ME); CS=os.path.dirname(SS)
sys.path.insert(0,SS); sys.path.insert(0,os.path.join(CS,"inside_visible"))
from ss_core import *
_c=os.getcwd(); from inside_core import FLOC; os.chdir(_c)
VAR=os.environ.get("AO","smooth"); GAP_US=0.5; GAP_RING=0.1
class AO(Model):
    def decode(self,p):
        p=p.reshape(self.n,4); u,lon,lat,s=p.T; isC=self.k=="C"
        M=np.zeros(self.n); Rv=np.zeros(self.n); Rc=np.zeros(self.n)
        M[isC]=10**(np.log10(MLO)+s[isC]*np.log10(MHI/MLO)); Rv[~isC]=RLO+s[~isC]*(RHI-RLO); Rc[~isC]=Rv[~isC]
        Rc[isC]=(3*M[isC]/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3)
        dmin=Rc+GAP_US; d=dmin+u*(RMAX-Rc-dmin)
        X=np.c_[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]*d[:,None]
        isV=self.k=="V"; M[isV]=-0.70*rho_gpc(X[isV])*vol(Rv[isV])*GPC**3
        return X,M,Rv,Rc
    def pattern(self,p):
        X,M,Rv,Rc=self.decode(p); return M@unit_field(X)-M.sum()*FLOC
    def violations(self,p):
        X,M,Rv,Rc=self.decode(p); d=np.linalg.norm(X,axis=1); isC=self.k=="C"
        v=[np.minimum(np.maximum(0,d+Rc-(RV-GAP_RING)),np.maximum(0,RV+GAP_RING-(d-Rc)))]       # no straddling of the ring
        if VAR=="smooth":
            rho=rho_gpc(X); Rs=np.where(isC,(3*M/(4*np.pi*rho*GPC**3))**(1/3),Rv)
            v.append(np.where(isC,np.maximum(0,np.log10(np.abs(M)/(MHI*rho/RHO))),0.0)*10)
            Rx=Rs
        else: Rx=Rc
        if self.n>1:
            i,j=np.triu_indices(self.n,1); v.append(np.maximum(0,Rx[i]+Rx[j]+GAP-np.linalg.norm(X[i]-X[j],axis=1)))
        return np.concatenate(v)
def encode(kinds,X,M,Rv):
    k=np.array(list(kinds)); X=np.array(X,float); M=np.array(M,float); Rv=np.array(Rv,float); p=np.zeros((len(k),4)); isC=k=="C"
    d=np.linalg.norm(X,axis=1); p[:,1]=np.arctan2(X[:,1],X[:,0]); p[:,2]=np.arcsin(np.clip(X[:,2]/d,-1,1))
    p[isC,3]=(np.log10(np.clip(M[isC],MLO,MHI))-np.log10(MLO))/np.log10(MHI/MLO); p[~isC,3]=(np.clip(Rv[~isC],RLO,RHI)-RLO)/(RHI-RLO)
    Rc=np.where(isC,(3*np.abs(M)/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3),Rv); dmin=Rc+GAP_US
    p[:,0]=np.clip((d-dmin)/(RMAX-Rc-dmin),0,1); return p.ravel()
