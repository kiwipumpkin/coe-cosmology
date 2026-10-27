# Phase5 routes (NEW Oct 26): Migkas A/A_all refit under the two background fixes. <H0>=68.64 fixed, age 13.787 Gyr (5% penalty as fitdk5),
# A=af^2/<af^2>, unbounded Nelder-Mead in log params; only k4b.sky physical/numerical guards (flagged if approached).
# Route A EPS: early-homogeneity constraint. LTB common bang + growing mode with the profile as the t0 density (dist2 physics, homog.drec);
#   require max_x |rho(t_rec,x)/rho(t_rec,1)-1| <= EPS over x in [0.2,1.8] (ray hit range). Hinge penalty 1e4*max(0,log10(dmax/EPS))^2.
#   Free: D, mc, l0, eta, L (k not imposed; implied local log-slope reported).
# Route B: local log-slope at our shell locked to CatWISE: -dlnrho/dlnx|_{x=1} = eta (D/L)^eta = 0.041037 -> L/D=(eta/0.041037)^(1/eta). Free: D, mc, l0, eta.
# Usage: python3 fitroute.py A 1e-3 | python3 fitroute.py B
import sys,os,json,numpy as np,warnings;warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/coepdf/joint_solve/phase4b');import k4b as K
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)));import homog as Hm
from scipy.optimize import minimize
O=os.path.dirname(os.path.abspath(__file__))+'/';R_=sys.argv[1];EPS=float(sys.argv[2]) if R_=='A' else None;TOL=0.05;SL=0.04103732750534994
TAG=f'A_{EPS:g}' if R_=='A' else 'B';XS=np.linspace(0.2,1.8,33);RLS=21.744
def unpack(x):
    D,mc,l0,eta=10**x[0],10**x[1],x[2],10**x[3]
    LD=10**x[4] if R_=='A' else (eta/SL)**(1/eta)
    return D,mc,l0,eta,LD
def parts(x):
    D,mc,l0,eta,LD=unpack(x);o=K.sky(D,1e-3,eta,mc,1.0,LD=LD)
    if o is None:return None
    a2=o[0]**2/np.mean(o[0]**2);m=K.cw.shifted(a2,l0%360)[K.cw.IMIG];chi=float(np.sum(K.cw.W*(m-K.cw.Y)**2))
    XR=np.linspace(max(0.02,1-RLS/D),1+RLS/D,33)  # (fix) x range actually reached by rays to last scattering (R_ls=21.74 Gpc fixed by H0+age)
    with np.errstate(all='ignore'):dr=Hm.drec(D,LD,eta,XR)[0]
    dmax=float(np.nanmax(np.abs(dr-1))) if np.all(np.isfinite(dr)) else np.inf
    return chi,o[1]['Ht']*K.TH,m,o,dmax
def obj(x):
    try:p=parts(x)
    except Exception:return 800.
    if p is None:return 800.
    f=p[0]+((p[1]/K.AGE-1)/TOL)**2
    if R_=='A':f+=1e4*max(0.,np.log10(max(p[4],1e-300)/EPS))**2 if np.isfinite(p[4]) else 1e6
    return f
if R_=='A':S=[[1.4368,-1.252,18.64,1.0607,np.log10(0.454)],[1.4368,-1.252,18.64,0.3,0.5],[1.44,-1.2,18.7,0.0,1.0],[1.44,-1.2,18.7,0.5,1.7]]
else:S=[[1.4368,-1.252,18.64,1.0607],[1.4459,-1.1866,18.81,0.3422],[1.443,-1.22,18.74,0.0],[1.44,-1.2,18.7,0.6]]
if os.environ.get('STARTS'):S=json.loads(os.environ['STARTS'])  # (added) override starts, e.g. warm-start eps=1e-3 from the eps=1e-2 optimum
ck=O+f'fit_{TAG}_ckpt.json';C=json.load(open(ck)) if os.path.exists(ck) else []
for s in S[len(C):]:
    r=minimize(obj,np.array(s),method='Nelder-Mead',options=dict(maxfev=int(os.environ.get('FEV','1500')),xatol=1e-5,fatol=1e-5,adaptive=True))
    C.append(dict(x0=s,x=r.x.tolist(),f=float(r.fun)));json.dump(C,open(ck,'w'));print(TAG,'start',round(r.fun,4),flush=True)
b=min(C,key=lambda c:c['f']);x=np.array(b['x'])
r=minimize(obj,x,method='Nelder-Mead',options=dict(maxfev=1500,xatol=1e-6,fatol=1e-6,adaptive=True));x=r.x if r.fun<b['f'] else x
chi,age,m,o,dmax=parts(x);Y=K.cw.Y;D,mc,l0,eta,LD=unpack(x)
R=dict(route=R_,eps=EPS,chi2=chi,obj=obj(x),r2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),D=D,mc=mc,l0=l0%360,eta=eta,LD=LD,L=LD*D,
 logslope_x1=float(eta*LD**-eta),catwise_logslope=SL,implied_k_local=float(eta*LD**-eta/SL),drec_max=dmax,age_Gyr=age,n=int(len(Y)),starts=[c['f'] for c in C],
 flags=[s for s,v in [('eta>500 (k4b guard 1e3)',eta>500),('mc<1e-5 (guard 1e-6)',mc<1e-5),('LD>50 (k4b guard 100)',LD>50),('LD<1e-5',LD<1e-5),('drec>EPS',R_=='A' and dmax>EPS*1.01)] if v])
json.dump(R,open(O+f'fit_{TAG}.json','w'),indent=1);print('DONE',R,flush=True)
