"""Shared pieces: Migkas data, soft-core density with CatWISE 3% lock (paper code), sky-circle observable."""
import json,os,sys,warnings,numpy as np
N='/workspace/cosmology-course/notes/'
_cwd=os.getcwd(); src=open(N+'catwise-density-levelset-migkas.py').read(); ns={'__name__':'x'}; exec(src[:src.index('def main')],ns); os.chdir(_cwd)
pc=ns['pc']; IMIG=ns['IMIG']; Y=ns['MA']; W=ns['W']; TH1=ns['TH1']; RV=ns['RV']; CHI2_NULL=ns['CHI2_NULL']
L_from_delta=ns['L_from_delta']; sigma_from_delta=ns['sigma_from_delta']; delta_of=ns['delta_of']; shifted=ns['shifted']
LDEG=ns['LDEG']; MLO,MHI=ns['MLO'],ns['MHI']; TH0=ns['TH0'] if 'TH0' in ns else None
import bubble_core as bc
TH0=bc.TH0; RG=pc.RG; CKMS=299792.458; C_GPC_GYR=0.306601
def r2u(y,m): return float(1-np.sum((y-m)**2)/np.sum((y-y.mean())**2))
def profile(DR,alpha):
    """Soft core rho ~ 1/(1+(r/L)^alpha) (beta=1, paper), L from CatWISE lock delta=0.031 at d=14.3 Gpc; truncated at Rb (paper rb_of 1.5).
    Returns enclosed-mass fraction table MT on RG (units of D, MT(1)=1), L/D, Rb."""
    D=DR*RV; sh=dict(alpha=alpha,beta=1.0)
    L=L_from_delta('soft',D,RV,sh); Rb=bc.rb_of(DR,1.5)
    S,MT=pc.tables('soft',dict(LD=L/D,alpha=alpha,beta=1.0),Rb)
    assert abs(delta_of('soft',D,L,RV,sh)-0.031)<1e-6
    return MT,L/D,Rb
def phi_table(MT):
    """Newtonian potential per unit (G M(<D)/D): -[m(r)/r + int_r^inf m'/r' dr'] (paper pc.phi_coe with mu=1)."""
    return lambda r: pc.phi_coe(np.atleast_1d(np.asarray(r,float)),1.0,MT)
def af_from_profile(vprof,DR,th,relativistic=False):
    rho=1.0/DR; t=np.radians(th); P=np.c_[1+rho*np.cos(t),rho*np.sin(t)]; r=np.hypot(*P.T)
    v=vprof(r)[:,None]*P/r[:,None]; vo=float(vprof(np.array([1.0]))[0])
    if not relativistic: vr=v-[vo,0]
    else:
        g=1/np.sqrt(1-vo**2); den=1-vo*v[:,0]; vr=np.c_[(v[:,0]-vo)/den,v[:,1]/(g*den)]
    return np.hypot(*vr.T)/rho
def score(af_on_grid,th,l0):
    """af sampled at th (deg, from anti-CoE). mean-normalise on 1-deg ring, shift by l0, compare with Migkas bins."""
    a=np.interp(TH1,np.append(th,360),np.append(af_on_grid,af_on_grid[0])); a=a/a.mean()
    m=shifted(a,l0)[IMIG]; return float(np.sum(W*(m-Y)**2)),m,shifted(a,l0)
