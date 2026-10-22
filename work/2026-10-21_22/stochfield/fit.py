import sys,json,numpy as np,healpy as hp
from scipy.optimize import minimize
from fwd import Fwd,T0
LMAX=int(sys.argv[1]);BG=sys.argv[2]=='bg'   # bg: include deterministic CoE background (smooth red + RS_cc) in ensemble mean
Pf=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt');Pf=Pf[(Pf[:,0]>=2)&(Pf[:,0]<=LMAX)]
lc=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt')[:,:2];l=np.arange(2,LMAX+1);Lt=np.interp(l,lc[:,0],lc[:,1])
z=np.load('/workspace/redweb2/realweb/los_flow_R/hl_R4_ns64.npz');ll=np.arange(LMAX+1)
Dfix=(ll*(ll+1)/(2*np.pi)*hp.anafast(hp.remove_dipole(z['sm']+z['rs_cc']),lmax=LMAX)*T0**2)[2:] if BG else np.zeros(len(l))
F=Fwd(lmax=LMAX)
def chi2(Dm,Dfield):
    e=np.where(Dm>Pf[:,1],Pf[:,3],Pf[:,2]);cv=2/(2*l+1)*Dfield**2
    return float(np.sum((Dm-Pf[:,1])**2/(e**2+cv))),float(np.sum((Dm-Pf[:,1])**2/e**2))
r2=lambda D:float(np.corrcoef(D,Pf[:,1])[0,1]**2)
NM=['lgA','n','lkc','chiw','lsw']
def P(x):return dict(lgA=x[0],n=x[1],kc=10**x[2],chiw=x[3],sw=10**x[4])
lk=np.log(F.k)
def reg(p):
    s=F.sig_delta(p);ramp=s**2                                   # L2 amplitude: rms delta (0.078 Gpc smoothing) squared
    w=F.k**3*F.Pphi(p);w=w/np.trapezoid(w,lk);d2=-4*(F.k/p['kc'])**2   # d2 lnP/dlnk2 of the parametric form
    rP=float(np.trapezoid(w*d2**2,lk))
    W,_=F.W(p);Wd=F.W(p,True);rR=float(F.R**2*np.sum(Wd**2)/max(np.sum(W**2),1e-300))   # radial roughness (R W'/W)^2
    return ramp+rP+rR,dict(amp=ramp,P=rP,rad=rR,sig_delta=s)
def model(x):
    p=P(x);Df=F.Dl(p);return Df+Dfix,Df
def obj(x,lam):
    if not(-30<x[0]<10 and -4<x[1]<6 and -1.5<x[2]<2.3 and -2<x[3]<10 and -1.5<x[4]<1.3):return 1e12
    Dm,Df=model(x);c=chi2(Dm,Df)[0]
    return c+lam*reg(P(x))[0] if lam>0 else c
def best(lam,starts):
    out=[]
    for x0 in starts:
        x0=np.array(x0,float);Dm,Df=model(x0);x0[0]+=np.log10(np.median(Pf[:,1])/max(np.median(Df),1e-300))   # level-start amplitude
        r=minimize(obj,x0,args=(lam,),method='Nelder-Mead',options=dict(maxiter=(1500 if LMAX<=30 else 800),xatol=1e-4,fatol=1e-4));out.append(r)
    return min(out,key=lambda r:r.fun)
rng=np.random.default_rng(1)
starts=[[0,n,0.3,cw,lsw] for n in (-2,1) for cw in (1,3.5) for lsw in (-0.3,0.3)]
res={'LMAX':LMAX,'BG':BG,'lcdm':dict(chi2_planck=chi2(Lt,0*Lt)[1],chi2_cv=chi2(Lt,Lt)[0],r2=r2(Lt)),'Lcurve':[]}
LAMS=[0]+list(np.geomspace(1e-2,1e3,11))
if len(sys.argv)>3:   # seed from a previous (l<=30) fit, reduced grid for cost
    prev=json.load(open(sys.argv[3]));starts=[e['x'] for e in prev['Lcurve'][::3]];LAMS=[0]+list(np.geomspace(1e-2,1e3,6))
r0=best(0,starts);x=r0.x
for lam in LAMS:
    r=best(lam,[x,r0.x]+([starts[1],starts[6]] if len(starts)>6 else []));Dm,Df=model(r.x);c=chi2(Dm,Df);g,gd=reg(P(r.x))
    res['Lcurve'].append(dict(lam=lam,x=r.x.tolist(),p=P(r.x),chi2_cv=c[0],chi2_planck=c[1],r2=r2(Dm),reg=g,**{'reg_'+k:v for k,v in gd.items()},D=Dm.tolist()))
    print('lam %.3g chi2cv %.2f chi2P %.2f r2 %.3f reg %.3g'%(lam,c[0],c[1],r2(Dm),g),{k:round(v,4) for k,v in P(r.x).items()},{k:'%.3g'%v for k,v in gd.items()},flush=True)
res['Dfix']=Dfix.tolist();res['planck']=Pf.tolist();res['lcdmD']=Lt.tolist()
json.dump(res,open(f'fit_l{LMAX}_{"bg" if BG else "nobg"}.json','w'),indent=1)
