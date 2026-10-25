# Replicates /workspace/cosmology-course/notes/exponential-migkas-fit.py (bubble_core.run geometry):
# sources on a ring at ONE observation distance R from us, in the plane containing the CoE;
# A/A_all(l) = a_f(l-l0)/<a_f>, a_f = H along direction at distance R; no latitude/region binning.
# Here a_f(psi;z) = c z / D_A(psi,z) from the LTB null-geodesic grid in ../zdist/runs (famZ.py). R free via z in [0.004,0.447].
import json,glob,os,csv,numpy as np,warnings;warnings.filterwarnings('ignore')
from scipy.optimize import minimize
from scipy.interpolate import RegularGridInterpolator as RGI
M=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
L=np.array([float(r['l_deg']) for r in M]);Y=np.array([float(r['A_over_Aall']) for r in M]);S=0.5*np.array([float(r['A_err_low'])+float(r['A_err_high']) for r in M])
zz=np.linspace(0.02,0.3,2000);WI=np.trapezoid(zz,zz)/np.trapezoid(zz**2,zz);C=299792.458
ZLO,ZHI=0.004,0.447;TH=np.arange(360.)
chi=lambda m:float(np.sum(((m-Y)/S)**2));R2=lambda m:float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2))
def setup(f):
    d=np.load(f);PS=d['PS'];Z=d['ZG'];DA=d['DA']
    P=np.r_[-PS[::-1],PS,360-PS[::-1]];Dx=np.r_[DA[::-1],DA,DA[::-1]]      # reflect psi about 0 and 180
    I=RGI((P,np.log(Z)),np.log(Dx),bounds_error=False,fill_value=None)
    def ring(z):  # a_f on 1-deg ring, psi=angle from away-CoE direction (bubble_core convention), mean-normalised
        ps=np.minimum(TH,360-TH);a=1/np.exp(I(np.c_[ps,np.full(360,np.log(z))]));return a/a.mean()
    return ring,d
def sh(a,l0): return np.interp((L-l0)%360,np.r_[TH,360],np.r_[a,a[0]])
def fit(f):
    ring,d=setup(f);o={}
    # G: exact replica (geometry only, linear): free z, l0
    best=None
    for z0 in np.geomspace(ZLO,ZHI,9):
        a=ring(z0);g=min(((chi(sh(a,l)),l) for l in range(0,360,4)))
        r=minimize(lambda p:1e9 if not ZLO<=10**p[0]<=ZHI else chi(sh(ring(10**p[0]),p[1])),[np.log10(z0),g[1]],method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-8))
        if best is None or r.fun<best.fun:best=r
    z,l0=10**best.x[0],best.x[1]%360;m=sh(ring(z),l0)
    o['G']=dict(chi2=chi(m),R2=R2(m),z=z,l0=l0,amp=float(ring(z).max()-ring(z).min()),A=m.tolist())
    # F: migkas_prof form with stretch term -> ring a_f(z): [a(l-l0)(1+bf WI cos)-vp/c WI sqrt(1-sb^2) cos(l-lpp)]^2
    LO=np.array([np.log10(ZLO),-360,-0.5,0,0,-1]);HI=np.array([np.log10(ZHI),720,0.5,3000,360,1])
    def mod(th):
        lz,l0,bf,vp,lpp,sb=th;t=np.radians(L-l0)
        return (sh(ring(10**lz),l0)*(1+bf*WI*np.cos(t))-vp/C*WI*np.sqrt(1-sb**2)*np.cos(np.radians(L-lpp)))**2
    pen=lambda th:1e9 if (np.any(th<LO) or np.any(th>HI)) else chi(mod(th))
    rng=np.random.default_rng(1);bf_=None
    st=[np.array([best.x[0],l0,0,1000,215,0]),np.array([best.x[0],l0,0,3000,215,0])]+[np.r_[LO[0]+(HI[0]-LO[0])*rng.random(),360*rng.random(),LO[2:]+(HI[2:]-LO[2:])*rng.random(4)] for _ in range(10)]
    for s in st:
        r=minimize(pen,s,method='Nelder-Mead',options=dict(maxiter=6000,xatol=1e-6,fatol=1e-8))
        if bf_ is None or r.fun<bf_.fun:bf_=r
    th=bf_.x.copy();m=mod(th);th[1]%=360
    o['F']=dict(chi2=chi(m),R2=R2(m),z=10**th[0],l0=th[1],bf=th[2],vp=th[3],lpp=th[4],sb=th[5],A=m.tolist())
    o['DR']=float(d['DR']);o['ratio']=float(d['ratio'])
    # R in Gpc: D_A at psi=90 times D, D = DR*14.30048 Gpc (famP convention, DRV)
    for k in 'GF':
        z=o[k]['z'];o[k]['DA_Gpc']=float(np.interp(np.log(z),np.log(d['ZG']),d['DA'][22])*d['DR']*14.30048287260589)
        for q in o[k]:
            if isinstance(o[k][q],np.floating):o[k][q]=float(o[k][q])
    return o
if __name__=='__main__':
    os.makedirs('res',exist_ok=True)
    for f in sorted(glob.glob('../zdist/runs/*.npz')):
        k=os.path.basename(f)[:-4];out='res/%s.json'%k
        if os.path.exists(out):continue
        o=fit(f);json.dump(o,open(out,'w'));print(k,round(o['G']['chi2'],3),round(o['G']['z'],4),round(o['G']['l0'],1),round(o['F']['chi2'],3),flush=True)
