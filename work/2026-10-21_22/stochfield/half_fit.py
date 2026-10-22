# Damped half-scale term: P = P0 k^n [1 + A cos(k rs+phi) + B e^{-k/kB} cos(k r2 + phi2)] e^{-(k/kD)^2}; r2 = rs/2 ('half') or free ('r2free')
# Checkpointed: each Nelder-Mead start writes half_ckpt_<mode>.json; rerun resumes.
import sys,os,json,numpy as np
from scipy.optimize import minimize
from multiprocessing import Pool
import acoustic_model as M
import acoustic_fit as A
def Pd(s,p,k=None):
    k=s.k if k is None else k;r2=p.get('r2',p['rs']/2)
    return 10**p['lgA']*k**p['n']*(1+p['Aw']*np.cos(k*p['rs']+p['phi'])+p.get('B',0)*np.exp(-k/p.get('kB',1e9))*np.cos(k*r2+p.get('phi2',0)))*np.exp(-(k/p['kD'])**2)
M.AFwd.Pd=Pd
J=json.load(open('acoustic_fit.json'));OLD=J['free_W']['best'];CORR0=np.array(J['free_W']['corr'])
NM=['lgA','n','Aw','rs','phi','lkD','chiw','lsw','B','phi2','lkB','r2']
def P(x,mode):
    p=dict(lgA=x[0],n=x[1],Aw=x[2],rs=x[3],phi=x[4]%(2*np.pi),kD=10**x[5],chiw=x[6],sw=10**x[7],B=x[8],phi2=x[9]%(2*np.pi),kB=10**x[10])
    if mode=='r2free':p['r2']=x[11]
    return p
def ok(x,mode):
    b=-25<x[0]<5 and -3<x[1]<6 and 0<=x[2]<=1 and 0.015<x[3]<0.25 and 0.5<x[5]<3.5 and 0.3<x[6]<=4.0 and np.log10(0.05)<x[7]<0.4 and 0<=x[8]<=3 and 0.5<x[10]<4.5
    return b and (mode!='r2free' or 0.01<x[11]<0.2)
CORR=CORR0.copy()
def obj(x,mode):
    if not ok(x,mode):return 1e15
    p=P(x,mode);Df=M.Dfield(p)*CORR;return A.chi2_fit(Df+M.Dfix)
def run(a):
    x0,mode,fev=a;r=minimize(obj,np.array(x0,float),args=(mode,),method='Nelder-Mead',options=dict(maxfev=fev,xatol=1e-5,fatol=1e-3,adaptive=True));return x0,r.x.tolist(),float(r.fun)
def full(x,mode):
    p=P(x,mode);Df=M.Dfield(p)*CORR;D=Df+M.Dfix
    return dict(x=list(x),p=p,chi2_fit=A.chi2_fit(D),**{k:A.stats(D,Df,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500),l30_2500=(30,2500)).items()})
if __name__=='__main__':
    mode=sys.argv[1];stage=sys.argv[2];ck=f'half_ckpt_{mode}.json'
    C=json.load(open(ck)) if os.path.exists(ck) else dict(starts=[],stage2=[],corr=None)
    if C['corr']:CORR[:]=C['corr']
    xo=OLD['x']
    if stage=='1':
        # seed B from residual amplitude: |Planck-CoE| in muK2 / D_CoE at first swings, phi2 grid incl. residual phase 5.55
        R=np.load('resid_free_W.npy');sig=np.interp(R[:,0],A.Pb[:,0],A.Pb[:,2]);Dold=M.Dfield(OLD['p'])*CORR0+M.Dfix
        m=R[:,0]<700;B0=float(np.median(np.abs(R[m,1])*sig[m]/np.interp(R[m,0],np.arange(len(Dold)),Dold)))*2
        print('B0 seed',B0,flush=True)
        seeds=[]
        for ph2 in (5.55,5.55-np.pi,0.0,np.pi/2,np.pi,3*np.pi/2):
            for B in (B0,0.5*B0):
                s=xo+[B,ph2,np.log10(250)]+([xo[3]/2] if mode=='r2free' else [])
                if mode=='r2free':s[11]*= (1.0 if B==B0 else 1.05)
                seeds.append(s)
        done=[tuple(np.round(s[0],6)) for s in C['starts']];todo=[s for s in seeds if tuple(np.round(s,6)) not in done]
        with Pool(8) as pl:
            for x0,x,f in pl.imap_unordered(run,[(s,mode,1500) for s in todo]):
                C['starts'].append([x0,x,f]);json.dump(C,open(ck,'w'));print('start',round(f,1),np.round(x,4).tolist(),flush=True)
    if stage=='2':   # refresh Limber->exact correction at current best, polish top 4
        C['starts'].sort(key=lambda t:t[2]);best=C['starts'][0][1]
        c,nn,rr=M.corr_curve(P(best,mode));CORR[:]=c;C['corr']=c.tolist();print('corr',np.round(rr,3).tolist(),flush=True)
        with Pool(4) as pl:
            for x0,x,f in pl.imap_unordered(run,[(s[1],mode,2500) for s in C['starts'][:4]]):
                C['stage2'].append([x0,x,f]);json.dump(C,open(ck,'w'));print('polish',round(f,1),flush=True)
    if stage=='3':
        b=min(C['stage2'] or C['starts'],key=lambda t:t[2])[1];C['best']=full(b,mode);json.dump(C,open(ck,'w'),indent=1);print(C['best'])
