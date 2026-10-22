# Resumable driver with checkpoints (acoustic_progress.json). Seeds free_W from the interrupted run's top starts.
import os;os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
import json,numpy as np
from multiprocessing import Pool
import acoustic_fit as A
from acoustic_model import Dfield,Dfix,EX,corr_curve
CK='acoustic_progress.json';S=json.load(open(CK)) if os.path.exists(CK) else {}
def save():json.dump(S,open(CK,'w'),indent=1)
def runc(a):
    x0,fix,lam,it,corr=a;A.CORR[:]=corr;return A.run((x0,fix,lam,it))
if 'lcdm' not in S:S['lcdm']={k:A.stats(A.Lt,A.Lt,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500),l30_2500=(30,2500)).items()};S['lcdm']['fit_measure_chi2']=A.chi2_fit(A.Lt);save()
one=np.ones(A.LMAX+1).tolist()
if 'fixed_shell_starts' not in S:
    seeds=[[-9,0.5,0.6,rs,ph,np.log10(400)] for rs in (0.042,0.063,0.084,0.105) for ph in (0,np.pi/2,np.pi,3*np.pi/2)]
    with Pool(8) as pl:R=pl.map(runc,[(s,True,0,0,one) for s in seeds])
    S['fixed_shell_starts']=sorted(R,key=lambda t:t[1]);save();print('fixed starts',S['fixed_shell_starts'][:3],flush=True)
S.setdefault('free_W_starts',[[[-8.363,-0.174,0.197,0.082,7.623,2.515,3.946,-1.301],10128.2],[[-8.37,-0.173,0.197,0.082,1.341,2.515,3.945,-1.289],10128.7],[[-6.592,-0.147,0.202,0.061,1.378,2.633,2.979,-1.301],10278.4],[[-6.433,-0.134,0.207,0.054,1.381,2.685,2.632,-1.267],10360.0]]);save()
for nm,fix in (('free_W',False),('fixed_shell',True)):
    T=S.setdefault(nm,{});R=S[nm+'_starts']
    while len(T.get('corr_hist',[]))<2:
        best=T['best_x'] if 'best_x' in T else R[0][0];c,nn,rr=corr_curve(A.P(best,fix))
        with Pool(4) as pl:RR=pl.map(runc,[(x,fix,0,1,c.tolist()) for x,_ in R[:4]])
        bx,f=min(RR,key=lambda t:t[1]);T['best_x']=bx;T.setdefault('corr_hist',[]).append(dict(nodes=nn.tolist(),ratio=rr.tolist(),chi2=f));save();print(nm,'corr',f,np.round(rr,3).tolist(),flush=True)
    if 'corr' not in T:c,_,_=corr_curve(A.P(T['best_x'],fix));T['corr']=c.tolist();save()
    c=np.array(T['corr']);T.setdefault('Lcurve',[]);done={e['lam'] for e in T['Lcurve']}
    todo=[l for l in (0,10,100,1000) if l not in done]
    if todo:
        with Pool(8) as pl:RR=pl.map(runc,[(x,fix,lam,1,T['corr']) for lam in todo for x in (T['best_x'],R[1][0])])
        for i,lam in enumerate(todo):
            x,f=min(RR[2*i:2*i+2],key=lambda t:t[1]);p=A.P(x,fix);Df=Dfield(p)*c;D=Df+Dfix
            T['Lcurve'].append(dict(lam=lam,x=x,p=p,obj=f,chi2_fit=A.chi2_fit(D),chi2_fit_cv=A.chi2_fit(D,Df,True),reg=A.reg(p),sig_delta=EX.sig_delta(p),**{k:A.stats(D,Df,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500),l30_2500=(30,2500)).items()}))
            print(nm,'lam',lam,round(T['Lcurve'][-1]['chi2_fit'],1),flush=True)
        save()
out=dict(params=A.NM,fixedW=A.FIXW,lcdm=S['lcdm'])
for nm in ('free_W','fixed_shell'):out[nm]=dict(starts=S[nm+'_starts'],corr_hist=S[nm]['corr_hist'],corr=S[nm]['corr'],Lcurve=sorted(S[nm]['Lcurve'],key=lambda e:e['lam']))
json.dump(out,open('acoustic_fit.json','w'),indent=1);print('DONE')
