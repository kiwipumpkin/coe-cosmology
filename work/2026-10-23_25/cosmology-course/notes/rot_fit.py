import numpy as np, pandas as pd, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.stats import chi2 as C, norm
N='/workspace/cosmology-course/notes/'
src=open(N+'catwise-density-levelset-migkas.py').read(); ns={'__name__':'x'}; exec(src[:src.index('Y = MA')],ns)
MA,MLO,MHI,IMIG=ns['MA'],ns['MLO'],ns['MHI'],ns['IMIG']; W=1/(0.5*(MLO+MHI))**2; Y=MA; C0=ns['CHI2_NULL']
Ld=np.loadtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',skiprows=1,usecols=0)
files={'exp':'catwise-density-levelset-migkas-exp_fixed15.2.csv','pl':'catwise-density-levelset-migkas-pl_fixed15.2.csv',
 'soft':'catwise-density-levelset-migkas-soft_fixed15.2.csv','eta':'catwise-eta-migkas-fit-d14.3.csv'}
name={'exp':r'$\rho_0\,e^{-r/L}$','pl':r'$\rho_0(1+r/L)^{-1.2}$','soft':r'$\rho_0/(1+(r/L)^{0.3})$','eta':r'$\rho_0\,e^{-(r/L)^{0.05}}$ (ℓ₀=19.9°)'}
r2=lambda m:1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2); c2=lambda m:float(np.sum(W*(m-Y)**2))
def sig(d,k): return float(norm.isf(C.sf(d,k))) if d>0 else 0.
out={};fig,ax=plt.subplots(figsize=(16,9))
ax.errorbar(Ld,MA,yerr=[MLO,MHI],fmt='o',color='k',ms=4,alpha=.6,label='Migkas bins')
for k,f in files.items():
    q=pd.read_csv(N+f); g=q.A_over_Aall.values; m0=g[IMIG]; lr=np.radians(q.ell_deg.values)
    rot=lambda p,l: 1-p[0]*np.sin(l+p[1])
    best=min((least_squares(lambda p:np.sqrt(W)*(m0*rot(p,lr[IMIG])-Y),[a,ph]) for a in(-.05,.05) for ph in np.linspace(0,6,7)),key=lambda r:r.cost)
    A,ph=best.x
    if A<0: A,ph=-A,ph+np.pi
    ph=ph%(2*np.pi); m1=m0*rot([A,ph],lr[IMIG])
    b,a=c2(m0),c2(m1)
    out[k]=dict(R2_before=r2(m0),R2_after=r2(m1),chi2_before=b,chi2_after=a,dchi2_rot=b-a,
      sigma_null_before=None,sigma_null_after_extra2=sig(C0-a,2),sigma_rot_only_2dof=sig(b-a,2),Amp=A,phi_deg=np.degrees(ph))
    ax.plot(q.ell_deg,g*rot([A,ph],lr),lw=3,label=f"{name[k]} × [1−A sin(ℓ+φ)]\nA={A:.3f}, φ={np.degrees(ph):.0f}°, R²={r2(m1):.3f}"); _x=(f"{A:.3f}, φ={np.degrees(ph):.0f}°, R²={r2(m1):.3f}")
ax.set_xlabel('Galactic longitude ℓ (deg)',fontsize=18); ax.set_ylabel('A / A_all',fontsize=18); ax.set_xlim(0,360)
ax.tick_params(labelsize=15); ax.grid(alpha=.3); ax.legend(fontsize=13,loc='upper left',bbox_to_anchor=(1.01,1))
ax.set_title('Migkas A/A_all, four CoE laws (fixed) + observer-rotation\nfactor 1−A sin(ℓ+φ), only A, φ fitted',fontsize=18,wrap=True)
fig.tight_layout(); fig.savefig(N+'catwise-AoverAall-four-laws-rotation.png',dpi=130,bbox_inches='tight')
print(C0,len(Y)); print(json.dumps(out,indent=1)); json.dump(out,open(N+'catwise-AoverAall-four-laws-rotation.json','w'),indent=1)
