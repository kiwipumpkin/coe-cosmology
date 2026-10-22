import sys,json,numpy as np
from scipy.optimize import minimize
from multiprocessing import Pool
from acoustic_model import *
Pf=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt');Pf=Pf[(Pf[:,0]>=2)&(Pf[:,0]<=LMAX)]
Pb=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-binned_R3.01.txt');lo=30+30*np.arange(len(Pb));hi=np.minimum(lo+29,LMAX)
th=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt');Lt=np.zeros(LMAX+1);Lt[2:]=np.interp(np.arange(2,LMAX+1),th[:,0],th[:,1])
low=Pf[Pf[:,0]<=29]
def binned(D):return np.array([D[a:b+1].mean() for a,b in zip(lo,hi)])
def chi2_fit(D,Df=None,cv=False):   # low l unbinned (2-29) + binned l>=30 ; Planck errors (side-dependent), optional cosmic variance of field
    l=low[:,0].astype(int);d=D[l];e=np.where(d>low[:,1],low[:,3],low[:,2])**2
    if cv:e=e+2/(2*l+1)*Df[l]**2
    c=np.sum((d-low[:,1])**2/e);db=binned(D);eb=Pb[:,2]**2
    if cv:eb=eb+2/((2*Pb[:,0]+1)*30)*binned(Df)**2
    return float(c+np.sum((db-Pb[:,1])**2/eb))
def stats(D,Df,a,b):   # unbinned full file, l in [a,b]
    m=(Pf[:,0]>=a)&(Pf[:,0]<=b);l=Pf[m,0].astype(int);d=D[l];e=np.where(d>Pf[m,1],Pf[m,3],Pf[m,2])**2
    c=float(np.sum((d-Pf[m,1])**2/e));ccv=float(np.sum((d-Pf[m,1])**2/(e+2/(2*l+1)*Df[l]**2)));return dict(chi2_planck=c,chi2_cv=ccv,r2=float(np.corrcoef(d,Pf[m,1])[0,1]**2),npts=int(m.sum()))
NM=['lgA','n','Aw','rs','phi','lkD','chiw','lsw']
FIXW=dict(chiw=4.0,sw=0.1)
def P(x,fix):
    p=dict(lgA=x[0],n=x[1],Aw=x[2],rs=x[3],phi=x[4]%(2*np.pi),kD=10**x[5])
    if fix:p.update(FIXW)
    else:p.update(chiw=x[6],sw=10**x[7])
    return p
def ok(x,fix):
    b=-25<x[0]<5 and -3<x[1]<6 and 0<=x[2]<=1 and 0.015<x[3]<0.25 and 0.5<x[5]<3.5
    return b and (fix or (0.3<x[6]<=4.0 and np.log10(0.05)<x[7]<0.4))
def reg(p):return p['Aw']**2+(0.1/p['sw'])**2+EX.sig_delta(p)**2
CORR=np.ones(LMAX+1)
def obj(x,fix,lam):
    if not ok(x,fix):return 1e15
    p=P(x,fix);Df=Dfield(p)*CORR;return chi2_fit(Df+Dfix)+lam*reg(p)
def run(a):
    x0,fix,lam,it=a;x0=np.array(x0,float)
    if lam==0 and it==0:   # level the amplitude first
        p=P(x0,fix);Df=Dfield(p);x0[0]+=np.log10(np.median(Pb[:,1])/max(np.median(binned(Df)),1e-300))
    r=minimize(obj,x0,args=(fix,lam),method='Nelder-Mead',options=dict(maxiter=2500,maxfev=4000,xatol=1e-5,fatol=1e-3,adaptive=True));return r.x.tolist(),float(r.fun)
if __name__=='__main__':
    out={'params':NM,'fixedW':FIXW,'lcdm':{k:stats(Lt,Lt,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500),l30_2500=(30,2500)).items()},'lcdm_fitmeasure':chi2_fit(Lt)}
    out['lcdm']['fit_measure_chi2']=out['lcdm_fitmeasure']
    print('LCDM',out['lcdm'],flush=True)
    res={}
    for fix in (False,True):
        nm='fixed_shell' if fix else 'free_W';seeds=[]
        for rs in (0.042,0.063,0.084,0.105):
            for ph in (0,np.pi/2,np.pi,3*np.pi/2):
                for W in ([None] if fix else [(4.0,np.log10(0.1)),(3.5,np.log10(0.3))]):
                    x=[-9,0.5,0.6,rs,ph,np.log10(400)]+([] if fix else list(W));seeds.append(x)
        CORR[:]=1
        with Pool(8) as pl:R=pl.map(run,[(s,fix,0,0) for s in seeds])
        R.sort(key=lambda t:t[1]);print(nm,'top starts',[(round(f,1),np.round(x,3).tolist()) for x,f in R[:4]],flush=True)
        # Limber -> exact-node correction iteration (2 passes)
        best=R[0][0];hist=[]
        for it in range(2):
            c,nn,rr=corr_curve(P(best,fix));CORR[:]=c
            with Pool(4) as pl:RR=pl.map(run,[(x,fix,0,1) for x,_ in R[:4]])
            best,f=min(RR,key=lambda t:t[1]);hist.append(dict(nodes=nn.tolist(),ratio=rr.tolist(),chi2=f));print(nm,'corr pass',it,f,np.round(rr,3).tolist(),flush=True)
        c,nn,rr=corr_curve(P(best,fix));CORR[:]=c
        Lc=[]
        for lam in (0,10,100,1000):
            with Pool(2) as pl:RR=pl.map(run,[(best,fix,lam,1),(R[1][0],fix,lam,1)])
            x,f=min(RR,key=lambda t:t[1]);p=P(x,fix);Df=Dfield(p)*CORR;D=Df+Dfix
            Lc.append(dict(lam=lam,x=x,p=p,obj=f,chi2_fit=chi2_fit(D),chi2_fit_cv=chi2_fit(D,Df,True),reg=reg(p),sig_delta=EX.sig_delta(p),
                **{k:stats(D,Df,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500),l30_2500=(30,2500)).items()}))
            print(nm,'lam',lam,round(Lc[-1]['chi2_fit'],1),{k:round(v,4) for k,v in p.items()},flush=True)
        res[nm]=dict(starts=[(x,f) for x,f in R],corr_hist=hist,corr=CORR.tolist(),Lcurve=Lc)
        json.dump(dict(out,**res),open('acoustic_fit.json','w'),indent=1)
