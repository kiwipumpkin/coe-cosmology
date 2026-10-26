# Phase4b (NEW): joint fit with no arbitrary d floor. <H0>=68.64 locked (age=Ht*977.8/68.64), age anchor 13.787 Gyr (TOL 5% penalty), k free.
# Bounds (log10 unless noted): D [2,316] Gpc, d [1e-3,60] Gpc, mc [1e-5,1e3], l0 [0,360] deg, eta [1e-4,1e3], k [1e-3,1e5].
# Physical guards only (k4b.sky): d<D/1.05, L root exists, vrmin>0. Checkpoints DE population each generation (de_ckpt.npz).
import numpy as np,json,os,sys,multiprocessing as mp
from scipy.optimize import differential_evolution,minimize
sys.path.insert(0,'/workspace/coepdf/joint_solve/phase4b');import k4b as K
TOL=0.05;O='/workspace/coepdf/joint_solve/phase4b/'
B=[(np.log10(2),2.5),(-3,np.log10(60)),(-5,3),(0,360),(-4,3),(-3,5)]
NAMES=['logD','logd','logmc','l0','logeta','logk']
def parts(x):
    D,d,mc,l0,eta,k=10**x[0],10**x[1],10**x[2],x[3],10**x[4],10**x[5];o=K.sky(D,d,eta,mc,k)
    if o is None:return None
    a2=o[0]**2/np.mean(o[0]**2);m=K.cw.shifted(a2,l0)[K.cw.IMIG];chi=float(np.sum(K.cw.W*(m-K.cw.Y)**2));age=o[1]['Ht']*K.TH
    return chi,age,m,o
def obj(x):
    import warnings;warnings.filterwarnings('ignore')
    if any(not(lo<=v<=hi) for v,(lo,hi) in zip(x,B)):return 800.
    try:p=parts(x)
    except Exception:return 800.
    if p is None:return 800.
    return p[0]+((p[1]/K.AGE-1)/TOL)**2
def rep(x):
    chi,age,m,o=parts(x);Y=K.cw.Y;D,d,mc,l0,eta,k=10**x[0],10**x[1],10**x[2],x[3],10**x[4],10**x[5]
    return dict(chi2=chi,obj=obj(x),r2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),D=D,d=d,L=o[2]*D,eta=eta,mc=mc,l0=l0,k=k,delta=0.031*k,Ht=o[1]['Ht'],age_Gyr=age,n=int(len(Y)))
if __name__=='__main__':
    ck=O+'de_ckpt.npz';NP=12*6
    rng=np.random.default_rng(7);init=np.array([[rng.uniform(lo,hi) for lo,hi in B] for _ in range(NP)])
    # seed a few members at phase4 optimum region with smaller d (does not constrain search)
    for i,ld in enumerate([0.301,0,-0.5,-1,-2,-3]):init[i]=[1.4459,ld,-1.1866,18.81,0.3422,1.9997]
    gen=0
    if os.path.exists(ck):z=np.load(ck);init=z['pop'];gen=int(z['gen']);print('resume gen',gen,flush=True)
    def cb(xk,convergence=None,_s=[gen]):
        _s[0]+=1;print('gen',_s[0],round(obj(xk),4),np.round(xk,3).tolist(),flush=True)
        json.dump(dict(gen=_s[0],x=list(map(float,xk)),f=float(obj(xk))),open(O+'de_best.json','w'))
    # population checkpoint via wrapper: save every generation from inside solver
    from scipy.optimize._differentialevolution import DifferentialEvolutionSolver as S
    with mp.get_context('spawn').Pool(8) as p:
        s=S(obj,B,popsize=12,maxiter=int(os.environ.get('MAXIT','120'))-gen,tol=1e-7,rng=4+gen,workers=p.map,updating='deferred',polish=False,init=init)
        for i,_ in enumerate(s):
            np.savez(ck,pop=s._scale_parameters(s.population),gen=gen+i+1)
            cb(s.x)
            if s.convergence<1e-7 or i+1>=s.maxiter:break
        x=s.x
    r2=minimize(obj,x,method='Nelder-Mead',options=dict(maxfev=3000,xatol=1e-6,fatol=1e-6,adaptive=True));x=r2.x if r2.fun<obj(x) else x
    R=rep(x);at=[NAMES[i] for i,(v,(lo,hi)) in enumerate(zip(x,B)) if NAMES[i]!='l0' and min(v-lo,hi-v)<0.02*(hi-lo)]
    json.dump(dict(bounds=dict(zip(NAMES,B)),best_x=list(map(float,x)),best=R,at_bound=at),open(O+'joint4b.json','w'),indent=1);print('DONE',R,at,flush=True)
