# NEW phase3 (audit E13): prof_core a_f=|v_rel|/r is H-like, but the locked fit compared a_f/<a_f> directly to Migkas A/A_all.
# Correct relation A/A_all=(H/<H>)^2 -> model A = a_f^2/<a_f^2>. (a) score the OLD locked curve as H-like; (b) refit {D,mc,l0,eta} with the squared model, L locked (d=14.3, delta=0.031).
import sys,os,json,importlib.util,numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
sp=importlib.util.spec_from_file_location('cw','/workspace/cosmology-course/notes/catwise-eta-density-migkas.py');cw=importlib.util.module_from_spec(sp);sp.loader.exec_module(cw)
DL=cw.RV
def pred(D,eta,mc,l0):
    o=cw.sky_locked(D,eta,mc,DL)
    if o is None:return None
    af=o[0];a2=af**2/np.mean(af**2);return cw.shifted(a2,l0),af,o[1]
def chi2(x):
    D,mc,l0,eta=10**x[0],10**x[1],x[2],10**x[3];p=pred(D,eta,mc,l0)
    return 800.0 if p is None else float(np.sum(cw.W*(p[0][cw.IMIG]-cw.Y)**2))
def run(x0):r=minimize(chi2,x0,method='Nelder-Mead',options=dict(maxfev=int(os.environ.get('FEV','600')),xatol=1e-4,fatol=1e-4));return list(map(float,r.x)),float(r.fun)
if __name__=='__main__':
    out={}
    # (a) old locked csv read as H/<H>
    c=np.loadtxt('/workspace/cosmology-course/notes/catwise-eta-migkas-fit-d14.3.csv',delimiter=',',skiprows=1);lc,af=c[:,0],c[:,1]
    A2=af**2/np.mean(af**2);m=np.interp(cw.TH1 if len(cw.TH1)==len(lc) else lc,lc,A2)
    p1=json.load(open('/workspace/coepdf/joint_solve/phase2/p1_migkas_D34.96.json'))['34.96']
    S=[[np.log10(34.96),np.log10(p1['mc']),p1['l0'],np.log10(max(p1['eta'],0.031))]]+[[np.log10(d),np.log10(m_),l,np.log10(e)] for d,m_,l,e in [(58.2,0.2,15,1),(46,4.5,21,1),(40,.5,14,.7),(30,1,15,1.5),(50,2,16,.5),(35,.8,18,.3),(28,.2,18.6,.8),(20,1,18,1)]]
    with Pool(int(os.environ.get('NP','4'))) as pl:R=pl.map(run,S)
    R.sort(key=lambda t:t[1]);x,f=R[0];D,mc,l0,eta=10**x[0],10**x[1],x[2],10**x[3];A,afb,rr=pred(D,eta,mc,l0)
    Y=cw.Y;mm=A[cw.IMIG];r2=float(1-np.sum((Y-mm)**2)/np.sum((Y-Y.mean())**2))
    out['refit_sq']=dict(D=D,mc=mc,l0=l0,eta=eta,chi2_A=f,r2_A=r2,r2_pearson=float(np.corrcoef(Y,mm)[0,1]**2),age_Gyr=float(rr['age']),Ht=float(rr['Ht']),
        H_over_mean_ptp=float(np.ptp(np.sqrt(A))/np.mean(np.sqrt(A))),starts=R,chi2_null=float(np.sum(cw.W*(Y-1)**2)))
    out['old_locked_chi2_A_linear']=9.058  # r2 json phase2/fig7 (A-space linear model, as originally fitted)
    np.savetxt('/workspace/coepdf/joint_solve/phase3/migkas_refit_sq_curve.csv',np.c_[cw.TH1,A,np.sqrt(A)/np.mean(np.sqrt(A))],delimiter=',',header='ell_deg,A_over_Aall_model,H_over_meanH_model',comments='')
    json.dump(out,open('/workspace/coepdf/joint_solve/phase3/migkas_refit_sq.json','w'),indent=1,default=float);print(json.dumps({k:v for k,v in out['refit_sq'].items() if k!='starts'},indent=1))
