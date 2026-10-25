import sys,os,json,numpy as np
from scipy.integrate import solve_ivp
BUB="/workspace/cosmology/papers/spreadsheet_model/bubble"; sys.path.insert(0,BUB); cwd=os.getcwd(); os.chdir(BUB)
import ptshell, prof_core as P
from bubble_core import IMIG, MA, MLO, MHI, TH1, shifted, geometry, rb_of, TH0, RKM, GPC
os.chdir(cwd)
RV=14.30048287260589; OS=0.049
w=json.load(open("migkas-geometry-fit.json"))["wls_fit"]; D,mcc,kin,kout,l0=[w[k] for k in ("D","mcc","kin","kout","l0")]
a=dict(DR=D/RV,Os=OS,mcc=mcc,rin=1-kin*RV/D,rout=1+kout*RV/D,Lc=0.005,mc=0.0,Ob=0.0)
rrN=P.run("ptshell",a,0.0,0.0); afN=rrN["af"]/rrN["af"].mean()
S,MT,m1=ptshell.build(a); DR=a["DR"]; q,eps=geometry(TH1,DR); taui=0.01
def radial(r0,v0,Mfun):
    s=solve_ivp(lambda t,y:[y[1],-Mfun(y[0])/y[0]**2],(taui,1.0),[r0,v0],rtol=1e-11,atol=1e-13)
    return s.y[0,-1],s.y[1,-1],s.y[0].min()
def runR(mode,Ht=0.95):
    r0=np.hypot(*q.T)*taui; v0=np.hypot(*q.T); u=q/np.hypot(*q.T)[:,None]
    for _ in range(3):
        mu=m1*Ht**2/2
        if mode=="newton":  # sanity: static Eulerian table, must reproduce P.run
            res=[radial(r,v,lambda x:mu*np.interp(x,P.RG,MT)) for r,v in zip(r0,v0)]
        else:  # LTB: each shell keeps the mass enclosed at its START radius (same r0, v0, M0 -> same E)
            res=[radial(r,v,(lambda M:lambda x:M)(mu*np.interp(r,P.RG,MT))) for r,v in zip(r0,v0)]
        R=np.array([x[0] for x in res]); V=np.array([x[1] for x in res]); Ht=V[0]/R[0]
    X=R[:,None]*u; Vv=V[:,None]*u; vrel=Vv[1:]-Vv[0]; af=np.hypot(*vrel.T)/eps
    return af/af.mean(),R,V,mu,Ht
afC,_,_,_,_=runR("newton"); afG,R,V,mu,Ht=runR("ltb")
print("newton self-check max diff",np.abs(afC-afN).max())
t0=Ht*TH0*3.156e16; Dm=D*GPC; ccode=2.998e5*t0/Dm
Mi=mu*np.interp(np.hypot(*q.T)*taui,P.RG,MT)
comp=2*Mi/R/ccode**2; compN=2*mu*np.interp(R,P.RG,MT)/R/ccode**2
print("c in code units",ccode,"Earth speed/c",V[0]/ccode,"max speed/c",V.max()/ccode)
print("max 2GM/Rc2 today LTB",comp.max(),"Earth",comp[0])
# along trajectory, initial
r0=np.hypot(*q.T)*taui; print("max 2GM/Rc2 at start (LTB masses)", (2*Mi/r0/ccode**2).max())
Y=MA; W=1/(0.5*(MLO+MHI))**2
fullN=shifted(afN,l0); fullG=shifted(afG,l0)
def st(m): mm=m[IMIG]; return float(np.sum(W*(mm-Y)**2)),float(1-np.sum((Y-mm)**2)/np.sum((Y-Y.mean())**2))
out=dict(chiN=st(fullN),chiG=st(fullG),maxdiff=float(np.abs(fullN-fullG).max()),Ht=Ht,age=Ht*TH0,ccode=ccode,maxcomp=float(comp.max()),earthcomp=float(comp[0]))
json.dump({k:(float(v) if not isinstance(v,tuple) else v) for k,v in out.items()},open("/tmp/gr_top_stats.json","w")); print(out); np.save("paper_repro/stretch/gr_top_curves.npy",np.vstack([fullN,fullG]))
print("GR curve min/max",fullG.min(),fullG.max(),"Newton",fullN.min(),fullN.max())
print("E per shell sign (LTB): min E", (0.5*(np.hypot(*q.T))**2 - Mi/(np.hypot(*q.T)*taui)).min())
print("R today range",R.min(),R.max(),"vs geometry",np.hypot(*q.T).min(),np.hypot(*q.T).max())
print("Mi range",Mi.min(),Mi.max(),"vs Newtonian M today range",(mu*np.interp(R,P.RG,MT)).min(),(mu*np.interp(R,P.RG,MT)).max())
print("max 2GM/Rc2 Newtonian table today", compN.max())
