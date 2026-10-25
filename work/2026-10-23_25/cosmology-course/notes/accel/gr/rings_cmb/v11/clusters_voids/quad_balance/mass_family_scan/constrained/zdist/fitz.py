# Fit Migkas 2020 Fig.6 A/A_all (71 bins, l=5..355, width dl=90, all b) using per-cluster LTB distances.
import sys,os,json,glob,csv,numpy as np,warnings;warnings.filterwarnings('ignore')
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.optimize import minimize
M=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
L=np.array([float(r['l_deg']) for r in M]);Y=np.array([float(r['A_over_Aall']) for r in M]);S=0.5*np.array([float(r['A_err_low'])+float(r['A_err_high']) for r in M])
K=list(csv.DictReader(open('/workspace/cosmology/migkas/migkas2020_tableC1.csv')))
zc=np.array([float(r['z']) for r in K]);lc=np.radians([float(r['l']) for r in K]);bc=np.radians([float(r['b']) for r in K])
nc=np.c_[np.cos(bc)*np.cos(lc),np.cos(bc)*np.sin(lc),np.sin(bc)]
INB=np.abs((np.degrees(lc)[None,:]-L[:,None]+180)%360-180)<=45   # 71 x 313 bin membership
zz=np.linspace(0.02,0.3,2000);WI=np.trapezoid(zz,zz)/np.trapezoid(zz**2,zz);C=299792.458
LO=np.array([0,-90,-0.5,0,0,-1]);HI=np.array([360,90,0.5,3000,360,1])
def setup(f):
    d=np.load(f);DA=d['DA'];Z=d['ZG'];PS=d['PS']
    w=np.sin(np.radians(PS));Dref=np.trapezoid(DA*w[:,None],PS,axis=0)/np.trapezoid(w,PS)
    I=RGI((PS,np.log(Z)),np.log(DA/Dref[None,:]),bounds_error=False,fill_value=None)
    return I,float(d['ratio'])
def aexp(I,lC,bC,Z=None):
    Z=zc if Z is None else Z
    lr,br=np.radians(lC),np.radians(bC);n=np.array([np.cos(br)*np.cos(lr),np.cos(br)*np.sin(lr),np.sin(br)])
    psi=180-np.degrees(np.arccos(np.clip(nc@n,-1,1)));psi=np.clip(psi,1,179)
    dl=-2*I(np.c_[psi,np.log(Z)])          # delta ln A_i = 2 ln(Dref/D_A)
    return np.exp(INB@dl/INB.sum(1)-dl.mean())
def model(th,I,Z=None):
    lC,bC,bf,vp,lpp,sb=th;t=np.radians(L-lC)
    return (np.sqrt(aexp(I,lC,bC,Z))*(1+bf*WI*np.cos(t))-vp/C*WI*np.sqrt(1-sb**2)*np.cos(np.radians(L-lpp)))**2
chi=lambda m:float(np.sum(((m-Y)/S)**2));R2=lambda m:float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2))
def fit(f,Z=None):
    I,rat=setup(f);out=dict(ratio=rat)
    # expansion only: grid then polish
    G=[(chi(aexp(I,a,b,Z)),a,b) for a in range(0,360,10) for b in range(-80,81,20)];G.sort()
    best=None
    for c0,a,b in G[:5]:
        r=minimize(lambda p:chi(aexp(I,p[0]%360,np.clip(p[1],-90,90),Z)),[a,b],method='Nelder-Mead',options=dict(xatol=1e-4,fatol=1e-7))
        if best is None or r.fun<best.fun:best=r
    lC,bC=best.x[0]%360,float(np.clip(best.x[1],-90,90));m=aexp(I,lC,bC)
    out['exp']=dict(chi2=chi(m),R2=R2(m),lC=lC,bC=bC,A=m.tolist(),amp=float(m.max()-m.min()))
    rng=np.random.default_rng(1);bestf=None
    pen=lambda th:1e9 if (np.any(th<LO) or np.any(th>HI)) else chi(model(th,I,Z))
    starts=[np.array([lC,bC,0,1000,215,0]),np.array([(lC+180)%360,-bC,0,3000,215,0])]+[LO+(HI-LO)*rng.random(6) for _ in range(14)]
    for s in starts:
        r=minimize(pen,s,method='Nelder-Mead',options=dict(maxiter=6000,xatol=1e-6,fatol=1e-8))
        if bestf is None or r.fun<bestf.fun:bestf=r
    m=model(bestf.x,I,Z);out['full']=dict(chi2=chi(m),R2=R2(m),theta=bestf.x.tolist(),A=m.tolist(),Aexp_part=aexp(I,bestf.x[0],bestf.x[1],Z).tolist())
    return out
if __name__=='__main__':
    R=json.load(open('fitz_results.json')) if os.path.exists('fitz_results.json') else {}
    for f in sorted(glob.glob('runs/*.npz')):
        k=os.path.basename(f)[:-4]
        if k in R: continue
        R[k]=fit(f);R[k]['z0']=fit(f,np.full_like(zc,0.003));print(k,R[k]['ratio'],R[k]['exp']['chi2'],R[k]['exp']['lC'],R[k]['exp']['bC'],R[k]['full']['chi2'],flush=True)
        json.dump(R,open('fitz_results.json','w'))
