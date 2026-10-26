# Step 2: TT spectrum. Replace the fixed CoE background D_fix (hl_R4 sm+rs_cc) by the spectrum of the shear-compensated sky
# from s2, then refit the 12-parameter acoustic + damped half-scale field (r2free) from the current best. Reports chi2 (Planck errors),
# Pearson r2, ensemble mean and seed-42 one-sky realisation, exactly as half_figs_fixed.py does.
import sys,os,json,numpy as np,healpy as hp
from scipy.optimize import minimize
OUT=os.path.abspath('.');os.chdir('/workspace/redweb2/stochfield');sys.path.insert(0,'.')
import half_fit as H,acoustic_fit as A,acoustic_model as M
T0=2.7255e6;comp=sys.argv[1];fev=int(sys.argv[2]) if len(sys.argv)>2 else 1500
C=json.load(open('fixed_ckpt_r2free.json'));H.CORR[:]=C['corr'];x0=np.array(C['best']['x']);ll=np.arange(M.LMAX+1)
def real(Df,seed=42):return Df*np.random.default_rng(seed).chisquare(2*np.maximum(ll,1)+1)/(2*np.maximum(ll,1)+1)
def report(Dfix,x):
    p=H.P(x,'r2free');Df=M.Dfield(p)*H.CORR;D=Df+Dfix;Dr=Dfix+real(Df)
    return dict(mean={k:A.stats(D,Df,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500)).items()},
                seed42={k:A.stats(Dr,Df,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500)).items()},chi2_fit=A.chi2_fit(D)),D,Dr,Df
old=M.Dfix.copy()
z=np.load(os.path.join(OUT,comp));Tot=z['Tot'];L=3*64-1
new=np.zeros(M.LMAX+1);new[:L+1]=ll[:L+1]*(ll[:L+1]+1)/(2*np.pi)*hp.anafast(hp.remove_dipole(Tot),lmax=L)*T0**2;new[:2]=0
r_old,*_=report(old,x0);r_new0,*_=report(new,x0)
def obj(x):
    if not H.ok(x,'r2free'):return 1e15
    p=H.P(x,'r2free');return A.chi2_fit(M.Dfield(p)*H.CORR+new)
res=minimize(obj,x0,method='Nelder-Mead',options=dict(maxfev=fev,xatol=1e-5,fatol=1e-3,adaptive=True))
r_new,D,Dr,Df=report(new,res.x)
o=dict(compensation=comp,baseline_oldDfix=r_old,newDfix_same_params=r_new0,newDfix_refit=r_new,x_refit=res.x.tolist(),p_refit=H.P(res.x,'r2free'),
       lcdm=J if (J:=H.J.get('lcdm')) else None,Dfix_old_l2_6=old[2:7].tolist(),Dfix_new_l2_6=new[2:7].tolist())
tag=os.path.basename(comp).replace('.npz','')
np.savez(os.path.join(OUT,f's3_{tag}.npz'),D=D,Dr=Dr,Df=Df,Dfix_new=new,Dfix_old=old,Pf=A.Pf,Lt=A.Lt)
json.dump(o,open(os.path.join(OUT,f's3_{tag}.json'),'w'),indent=1,default=float)
for k in ('baseline_oldDfix','newDfix_same_params','newDfix_refit'):print(k,json.dumps({a:{b:round(o[k][a][b]['chi2_planck'],1) for b in o[k][a]} for a in ('mean','seed42')}),'r2',round(o[k]['mean']['l2_2500']['r2'],4))
