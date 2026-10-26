# Phase4 (NEW): <H0> locked (H sky shape is H0-independent; H0 enters only age = Ht*977.8/68.64 Gyr).
# Model: prof_core 'eta' bubble (the Migkas sky + age code). A/A_all model = af^2/<af^2> (E13).
# CatWISE gradient lock delta = 0.031*k (k<1: CatWISE overestimates by 1/k; k>1: underestimates). 'free': L free, no lock.
# Free: D, d (CatWISE depth), mc, l0, eta (+L if free). Age modes: none / ±20% / ±5% of 13.787 Gyr (penalty).
import importlib.util,numpy as np,json,os,sys,itertools
from multiprocessing import Pool
from scipy.optimize import minimize
sp=importlib.util.spec_from_file_location('cw','/workspace/cosmology-course/notes/catwise-eta-density-migkas.py');cw=importlib.util.module_from_spec(sp);sp.loader.exec_module(cw)
H0=68.64;TH=977.8/H0;AGE=13.787
def sky(D,d,eta,mc,k,LD=None):
    if not(D>d*1.05 and 1e-4<mc<80 and 1e-3<eta<8 and 1<d<60):return None
    if LD is None:
        try:L=cw.L_from_delta_eta(D,d,eta,0.031*k)
        except Exception:return None
        LD=L/D
    if not(1e-4<LD<3):return None
    try:rr=cw.pc.run('eta',dict(DR=D/cw.RV,LD=LD,eta=eta),mc,0.0)
    except Exception:return None
    af=rr['af']
    if not np.all(np.isfinite(af)) or rr['vrmin']<=0:return None
    return af,rr,LD
def unpack(x,free):
    D,d,mc,l0,eta=10**x[0],10**x[1],10**x[2],x[3],10**x[4];LD=10**x[5] if free else None;return D,d,mc,l0,eta,LD
def obj(x,k,free,tol):
    D,d,mc,l0,eta,LD=unpack(x,free);o=sky(D,d,eta,mc,k,LD)
    if o is None:return 800.
    a2=o[0]**2/np.mean(o[0]**2);c=float(np.sum(cw.W*(cw.shifted(a2,l0)[cw.IMIG]-cw.Y)**2))
    if tol:age=o[1]['Ht']*TH;dev=abs(age/AGE-1)-tol;c+=1e4*max(dev,0)**2
    return c
def job(a):
    k,free,tol,x0=a;r=minimize(obj,x0,args=(k,free,tol),method='Nelder-Mead',options=dict(maxfev=int(os.environ.get('FEV','700')),xatol=1e-4,fatol=1e-4,adaptive=True))
    return k,free,tol,list(map(float,r.x)),float(r.fun)
if __name__=='__main__':
    S=[(34.75,14.3,3.9,20,0.031,.004),(25,14.3,0.2,20,1,.1),(50,10,0.2,15,0.3,.03),(20,8,1,20,2,.3),(40,20,0.5,25,0.1,.01),(15,6,0.1,18,0.5,1)]
    cases=[(k,False) for k in (0.1,0.2,0.5,1,2,5,10)]+[(1,True)]
    jobs=[]
    for (k,free),tol in itertools.product(cases,(0,0.2,0.05)):
        for s in S:
            x=[np.log10(s[0]),np.log10(s[1]),np.log10(s[2]),s[3],np.log10(s[4])]+([np.log10(s[5])] if free else []);jobs.append((k,free,tol,x))
    R={}
    with Pool(int(os.environ.get('NP','8'))) as p:
        for k,free,tol,x,f in p.imap_unordered(job,jobs):
            key=f"{'free' if free else k}|{tol}"
            if key not in R or f<R[key][1]:R[key]=(x,f)
            print(key,round(f,3),flush=True)
    out={}
    for key,(x,f) in R.items():
        kk,tol=key.split('|');free=kk=='free';k=1. if free else float(kk);tol=float(tol)
        D,d,mc,l0,eta,LD=unpack(x,free);o=sky(D,d,eta,mc,k,LD)
        if o is None:out[key]=dict(fail=True,f=f);continue
        a2=o[0]**2/np.mean(o[0]**2);m=cw.shifted(a2,l0)[cw.IMIG];Y=cw.Y
        chi=float(np.sum(cw.W*(m-Y)**2));H=np.sqrt(a2)
        out[key]=dict(k=None if free else k,free=free,age_tol=tol,chi2=chi,chi2_flat=float(np.sum(cw.W*(1-Y)**2)),r2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),
            D=D,d=d,L=LD*D,LD=LD,eta=eta,mc=mc,l0=l0,Ht=o[1]['Ht'],age_Gyr=o[1]['Ht']*TH,H_ptp=float(H.max()-H.min()),n=int(len(Y)))
    json.dump(out,open('/workspace/coepdf/joint_solve/phase4/k_scan.json','w'),indent=1)
