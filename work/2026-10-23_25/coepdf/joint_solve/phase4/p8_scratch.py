# Phase4 (NEW): stochastic field fitted FROM SCRATCH at the phase4 background. Env: TAG R DG DFIXF NDRAW NST NP FEV.
# Phase3 (NEW): p8 at R=R_ls(H0) from env R, Dfix from env DFIXF; start = phase2 p8 best rescaled by s_new/s_old. Out: phase3 p8_<TAG>.*
# Phase2 step 1: rebuild the TT stochastic field at the chosen consistent background D=34.96 Gpc, R=R_ls=eta_rec*D=0.18957 Gpc
# (fwd.kernels had D=42 hard-coded in the moving-lens factor 42*Rd and R=4). All lengths in P(k)/W scaled by s=R/4, wavenumbers by 1/s,
# as a starting point; then multistart Nelder-Mead (15 params, as p4), Dfix = low-l spectrum of the p5 run-a residual (bg+bodies at R).
# Env: FEV NST NP STAGE(fit|report)
import sys,os,json,numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
O='/workspace/coepdf/joint_solve/phase4/';P2='/workspace/coepdf/joint_solve/phase2/';TAG=os.environ['TAG'];os.chdir('/workspace/redweb2/stochfield');sys.path.insert(0,'.')
import fwd as F
DG=float(os.environ['DG']);R=float(os.environ['R']);s=R/4.0;K0=100.0/s
def kernelsD(R,dch):
    ch=np.arange(0,R,dch)+dch/2;te=F.t_of_chi(np.r_[ch-dch/2,R],R);tt=F.t_of_chi(ch,R);dt=-np.diff(te)
    ai=np.interp(np.log(tt),F.lTs,F.a_);Am=0.5*(F.A(te[:-1])+F.A(te[1:]))
    KI=(2/F.c**2)*(F.A(te[:-1])-F.A(te[1:]))/dch;Rd=F.Stf(np.c_[np.log(tt),np.ones_like(tt)])
    KM=(2/F.c**2)*(-(Am/ai*dt)*DG*Rd)/dch;SW=F.A(F.t_of_chi(R-1e-9,R))/(3*F.c**2);return ch,KI,KM,SW
F.kernels=kernelsD
import acoustic_model as M,acoustic_fit as A
def Pd(S,p,k=None):
    k=S.k if k is None else k;lk=np.log(k/K0)
    return 10**p['lgA']*k**(p['n']+p['alpha']*lk)*(1+p['Aw']*np.exp(p['beta']*lk)*np.cos(k*p['rs']+p['phi'])+p['B']*np.exp(-k/p['kB'])*np.cos(k*p['r2']+p['phi2']))*np.exp(-(k/p['kD'])**p['m'])
M.AFwd.Pd=Pd
M.EX=M.AFwd(R=R,lmax=M.LX,kmax=200/s,kmin=1e-2/s,nk=500,dch=0.004*s)
Dfix=np.zeros(M.LMAX+1);a=np.load(os.environ['DFIXF'])['Dl_resid'];Dfix[2:41]=a[2:41];M.Dfix=Dfix
CORR=np.ones(M.LMAX+1)
def P(x):return dict(lgA=x[0],n=x[1],Aw=x[2],rs=x[3],phi=x[4]%(2*np.pi),kD=10**x[5],chiw=x[6],sw=10**x[7],B=x[8],phi2=x[9]%(2*np.pi),kB=10**x[10],r2=x[11],alpha=x[12],beta=x[13],m=x[14])
ls=np.log10(s)
def ok(x):return (-40<x[0]<20 and -3<x[1]<6 and 0<=x[2]<=1 and 0.015*s<x[3]<0.25*s and 0.5-ls<x[5]<3.5-ls and 0.3*s<x[6]<=R and np.log10(0.05*s)<x[7]<0.4+ls
   and 0<=x[8]<=3 and 0.5-ls<x[10]<4.5-ls and 0.01*s<x[11]<0.2*s and -1<x[12]<1 and -3<x[13]<3 and 0.5<x[14]<6)
def obj(x):
    if not ok(x):return 1e15
    return A.chi2_fit(M.Dfield(P(x))*CORR+Dfix)
def run(x0):
    FEV=int(os.environ.get('FEV','3000'));r=minimize(obj,x0,method='Nelder-Mead',options=dict(maxfev=FEV,xatol=1e-6,fatol=1e-4,adaptive=True));return list(x0),r.x.tolist(),float(r.fun)
ll=np.arange(M.LMAX+1);g42=np.random.default_rng(42).chisquare(2*np.maximum(ll,1)+1)/(2*np.maximum(ll,1)+1)
def report(x):
    Df=M.Dfield(P(x))*CORR;D=Df+Dfix;Dr=Dfix+Df*g42;rg=dict(l2_30=(2,30),l30_2500=(30,2500),l2_2500=(2,2500))
    return dict(fit_measure=A.chi2_fit(D),mean={k:A.stats(D,Df,*r) for k,r in rg.items()},seed42={k:A.stats(Dr,Df,*r) for k,r in rg.items()}),D,Dr,Df
# FROM SCRATCH: no reuse of earlier field parameters. Uniform random draws inside ok() box, screened.
LO=np.array([-40,-3,0,0.015*s,0,0.5-ls,0.3*s,np.log10(0.05*s),0,0,0.5-ls,0.01*s,-1,-3,0.5]);HI=np.array([20,6,1,0.25*s,2*np.pi,3.5-ls,R,0.4+ls,3,2*np.pi,4.5-ls,0.2*s,1,3,6])
def draw(rng):
    while True:
        t=LO+(HI-LO)*rng.random(15)
        if ok(t):return t
def scr(seed):
    rng=np.random.default_rng(seed);T=[draw(rng) for _ in range(int(os.environ.get('NDRAW','150')))]
    # amplitude is degenerate-ish: optimise lgA on a 1-D grid for each draw
    out=[]
    for t in T:
        g=np.linspace(-40,20,31);f=[obj(np.r_[v,t[1:]]) for v in g];t[0]=g[int(np.argmin(f))];out.append((float(min(f)),t.tolist()))
    return out
xs=None
if __name__=='__main__':
    ck=O+f'p8_{TAG}_ckpt.json';C=json.load(open(ck)) if os.path.exists(ck) else dict(scan=None,starts=[],polish=[],corr=None)
    if C['corr']:CORR[:]=C['corr']
    st=os.environ.get('STAGE','fit')
    if C['scan'] is None:
        with Pool(int(os.environ.get('NP','8'))) as pl:Rs=sum(pl.map(scr,range(8)),[])
        Rs.sort(key=lambda t:t[0]);C['scan']=dict(top=Rs[:40]);json.dump(C,open(ck,'w'));print('screen best',[round(r[0],1) for r in Rs[:8]],flush=True)
    if st=='fit':
        NST=int(os.environ.get('NST','16'));starts=[np.array(t[1]) for t in C['scan']['top'][:NST]]
        with Pool(int(os.environ.get('NP','6'))) as pl:
            for a0,x,f in pl.imap_unordered(run,starts[len(C['starts']):]):
                C['starts'].append([a0,x,f]);json.dump(C,open(ck,'w'));print('start',round(f,2),flush=True)
        C['starts'].sort(key=lambda t:t[2]);c,nn,rr=M.corr_curve(P(np.array(C['starts'][0][1])));CORR[:]=c;C['corr']=c.tolist();print('corr',np.round(rr,3).tolist(),flush=True)
        with Pool(3) as pl:C['polish']=pl.map(run,[np.array(t[1]) for t in C['starts'][:3]])
        json.dump(C,open(ck,'w'))
    b=min(C['polish'] or C['starts'],key=lambda t:t[2])[1];rep,D,Dr,Df=report(np.array(b))
    json.dump(dict(D_Gpc=DG,R_Gpc=R,scale=s,best_x=b,best_p=P(np.array(b)),best=rep,scan=C['scan'],n_starts=len(C['starts'])),open(O+f'p8_{TAG}.json','w'),indent=1,default=float)
    np.savez(O+f'p8_{TAG}.npz',D=D,Dr=Dr,Df=Df,Dfix=Dfix,Lt=A.Lt,Pf=A.Pf)
    print('BEST mean',{k:round(rep['mean'][k]['chi2_planck'],2) for k in rep['mean']},'seed42',{k:round(rep['seed42'][k]['chi2_planck'],2) for k in rep['seed42']},flush=True)
