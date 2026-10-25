import sys,os,json,time,numpy as np
BUB="/workspace/cosmology/papers/spreadsheet_model/bubble"; sys.path.insert(0,BUB); _c=os.getcwd(); os.chdir(BUB)
import ptshell, prof_core as P
from bubble_core import IMIG, MA, MLO, MHI, TH1, shifted, geometry
sys.path.insert(0,"/workspace/cosmology/papers/spreadsheet_model/integrated")
from integ_core import best_sine, A_from
os.chdir(_c)
RV=14.30048287260589; OS=0.049; W=1/(0.5*(MLO+MHI))**2; Y=MA
def nsky(D,mcc,kin,kout):
    """Newtonian (Chris's integrator) = matched-today GR"""
    rin=1-kin*RV/D; rout=1+kout*RV/D
    if not (D>0 and mcc>0 and kin>1.001 and kout>1.001 and 0.008<rin<0.99): return None
    rr=P.run("ptshell",dict(DR=D/RV,Os=OS,mcc=mcc,rin=rin,rout=rout,Lc=0.005,mc=0.0,Ob=0.0),0.0,0.0)
    if not np.all(np.isfinite(rr["af"])) or rr["vrmin"]<=0: return None
    return rr["af"]/rr["af"].mean()
def age_exp(R,M,E):
    out=np.empty_like(R); p=E>0; k=~p
    A=M[p]/(2*E[p]); et=np.arccosh(1+R[p]/A); out[p]=M[p]/(2*E[p])**1.5*(np.sinh(et)-et)
    A=M[k]/(-2*E[k]); c=np.clip(1-R[k]/A,-1,1); et=np.arccos(c); out[k]=M[k]/(-2*E[k])**1.5*(et-np.sin(et))
    return out
def ltb(D,mcc,kin,kout,full=False):
    """Nested LTB dust shells. Today's density = core+shell table; M(R)=mass enclosed today, conserved by each shell;
    common Big Bang at t=0 (the code's Hubble-flow start); E(R) solved so each shell is at R at t=1 (today)."""
    if not (D>0 and mcc>0 and kin>1.001 and kout>1.001): return None
    rin=1-kin*RV/D; rout=1+kout*RV/D
    if not (0.008<rin<0.99): return None
    a=dict(DR=D/RV,Os=OS,mcc=mcc,rin=rin,rout=rout,Lc=0.005)
    q,eps=geometry(TH1,D/RV); S,MT,m1=ptshell.build(a); R=np.hypot(*q.T); u=q/R[:,None]; Ht=0.95
    for _ in range(4):
        mu=m1*Ht**2/2; M=np.maximum(mu*np.interp(R,P.RG,MT),1e-30)
        lo=-M/R*(1-1e-12); hi=np.full_like(R,1e4)
        if np.any(age_exp(R,M,lo)<1): return None
        for _ in range(100):
            m=0.5*(lo+hi); m=np.where(m==0,1e-300,m); g=age_exp(R,M,m)>1; lo=np.where(g,m,lo); hi=np.where(g,hi,m)
        E=0.5*(lo+hi); v=np.sqrt(np.maximum(2*M/R+2*E,0)); Ht=v[0]/R[0]
    V=v[:,None]*u; vrel=V[1:]-V[0]; af=np.hypot(*vrel.T)/eps
    if full: return dict(af=af/af.mean(),R=R,v=v,M=M,E=E,Ht=Ht,mu=mu)
    return af/af.mean()
def chi(m): mm=m[IMIG]; return float(np.sum(W*(mm-Y)**2))
def r2(m): mm=m[IMIG]; return float(1-np.sum((Y-mm)**2)/np.sum((Y-Y.mean())**2))
def sine(c,l0):
    th=TH1-l0; a,ph=best_sine(c[IMIG],th[IMIG],Y); return A_from(c,th,a,ph),a,ph
