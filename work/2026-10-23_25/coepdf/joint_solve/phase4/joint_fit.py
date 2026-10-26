# Phase4 (NEW): joint fit, <H0>=68.64 locked, gradient free & continuous.
# Params: log D, log d, log mc, l0, log eta, log k  where CatWISE-lock delta = 0.031*k (k = true/CatWISE gradient ratio).
# Objective: Migkas chi2 with A=af^2/<af^2> + age penalty ((age/13.787-1)/TOL)^2, age=Ht*977.8/68.64 Gyr (prof_core = the one dynamics code used).
# Steps: differential evolution (global) -> profile chi2(k) -> numerical Hessian correlations. Env TOL (default 0.05).
import numpy as np,json,os,sys
import multiprocessing as mp
from scipy.optimize import differential_evolution,minimize
sys.path.insert(0,'/workspace/coepdf/joint_solve/phase4');import k_scan as K
TOL=float(os.environ.get('TOL','0.05'));O=f'/workspace/coepdf/joint_solve/phase4/joint_T{TOL:g}'
B0=[(np.log10(10),np.log10(80)),(np.log10(2),np.log10(40)),(-3,np.log10(30)),(0,360),(-3,np.log10(8)),(-2,2)]
B=B0
def parts(x):
    D,d,mc,l0,eta,k=10**x[0],10**x[1],10**x[2],x[3],10**x[4],10**x[5];o=K.sky(D,d,eta,mc,k)
    if o is None:return None
    a2=o[0]**2/np.mean(o[0]**2);m=K.cw.shifted(a2,l0)[K.cw.IMIG];chi=float(np.sum(K.cw.W*(m-K.cw.Y)**2));age=o[1]['Ht']*K.TH
    return chi,age,m,o
def obj(x):
    if any(not(lo<=v<=hi) for v,(lo,hi) in zip(x,B)):return 800.
    p=parts(x)
    if p is None:return 800.
    return p[0]+((p[1]/K.AGE-1)/TOL)**2
NAMES=['logD','logd','logmc','l0','logeta','logk']
def prof(a):
    lk,x0=a
    f=lambda y:obj(np.r_[y,lk]);best=None
    for s in (x0,np.r_[x0[:2],x0[2]+.5,x0[3:]],np.r_[x0[0]+.15,x0[1:]],np.r_[x0[:4],x0[4]+.5]):
        r=minimize(f,s,method='Nelder-Mead',options=dict(maxfev=900,xatol=1e-4,fatol=1e-4,adaptive=True))
        if best is None or r.fun<best.fun:best=r
    return lk,list(map(float,best.x)),float(best.fun)
def rep(x):
    chi,age,m,o=parts(x);Y=K.cw.Y;D,d,mc,l0,eta,k=10**x[0],10**x[1],10**x[2],x[3],10**x[4],10**x[5]
    return dict(chi2=chi,obj=obj(x),r2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),D=D,d=d,L=o[2]*D,eta=eta,mc=mc,l0=l0,k=k,delta=0.031*k,Ht=o[1]['Ht'],age_Gyr=age,n=int(len(Y)))
if __name__=='__main__':
    with mp.get_context('spawn').Pool(8) as p:
        r=differential_evolution(obj,B,popsize=16,maxiter=int(os.environ.get('MAXIT','80')),tol=1e-6,seed=4,workers=p.map,updating='deferred',polish=False)
    x=r.x;r2=minimize(obj,x,method='Nelder-Mead',options=dict(maxfev=3000,xatol=1e-5,fatol=1e-5,adaptive=True));x=r2.x if r2.fun<r.fun else x
    out=dict(TOL=TOL,de_fun=float(r.fun),best_x=list(map(float,x)),best=rep(x));print(out['best'],flush=True)
    json.dump(out,open(O+'.json','w'),indent=1)
    G=np.linspace(-2,2,33)
    with mp.get_context('spawn').Pool(8) as p:P=sorted(p.map(prof,[(lk,x[:5]) for lk in G]))
    out['profile']=[dict(logk=lk,obj=f,**{k:v for k,v in rep(np.r_[y,lk]).items()}) for lk,y,f in P];json.dump(out,open(O+'.json','w'),indent=1)
    # Hessian (central differences) in the 6 fit coordinates
    h=np.array([.01,.01,.02,.5,.02,.02]);n=6;Hm=np.zeros((n,n));f0=obj(x)
    for i in range(n):
        for j in range(i,n):
            def F(a,b):e=np.zeros(n);e[i]+=a*h[i];e[j]+=b*h[j];return obj(x+e)
            Hm[i,j]=Hm[j,i]=(F(1,1)-F(1,-1)-F(-1,1)+F(-1,-1))/(4*h[i]*h[j])
    try:C=np.linalg.inv(Hm/2);sd=np.sqrt(np.abs(np.diag(C)));out['corr']=(C/np.outer(sd,sd)).tolist();out['sigma']=dict(zip(NAMES,sd.tolist()));out['hess_eig']=np.linalg.eigvalsh(Hm).tolist()
    except Exception as e:out['corr_err']=str(e)
    json.dump(out,open(O+'.json','w'),indent=1);print('done',flush=True)
