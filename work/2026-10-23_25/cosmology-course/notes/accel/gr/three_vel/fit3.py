import sys,os,numpy as np,json
from scipy.optimize import least_squares
B="/workspace/cosmology/papers/spreadsheet_model/bubble";sys.path.insert(0,B);c=os.getcwd();os.chdir(B)
from bubble_core import IMIG,MA,MLO,MHI;os.chdir(c)
L=IMIG.astype(float);Y=MA;S=0.5*(MLO+MHI)
zz=np.linspace(0.02,0.3,2000);winv=np.trapezoid(zz,zz)/np.trapezoid(zz**2,zz)  # <1/z>, uniform in volume
bp=0.00123;lp,bpl=264.021,48.253;cb=np.cos(np.radians(bpl))
def model(p,lpdir=lp,fixf=None):
    q,df,l0=p[:3]; df=df if fixf is None else fixf
    t=np.radians(L-l0); h=1+q*np.cos(t)**2+df*np.cos(t)-bp*winv*cb*np.cos(np.radians(L-lpdir))
    th=np.radians(np.arange(360)-l0);hb=np.mean(1+q*np.cos(th)**2)
    return (h/hb)**2
def fit(fixf=None):
    r=least_squares(lambda p:(model(p,fixf=fixf)-Y)/S,[0.08,0.0,17.8]);return r.x
st={}
for name,ff in [("full",None),("frame0",0.0)]:
    p=fit(ff);m=model(p,fixf=ff);chi=float(np.sum(((m-Y)/S)**2));R2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2))
    st[name]=dict(q=p[0],df=(p[1] if ff is None else 0),l0=p[2]%360,chi2=chi,R2=R2)
# errors full
r=least_squares(lambda p:(model(p)-Y)/S,[st['full'][k] for k in('q','df','l0')]);J=r.jac;cov=np.linalg.inv(J.T@J);st['err']=list(np.sqrt(np.diag(cov)));st['corr_q_df']=cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
st['flat_chi2']=float(np.sum(((Y-np.average(Y,weights=1/S**2))/S)**2));st['inv_z_eff']=winv;st['pec_dipole_H']=bp*winv*cb
f=st['full'];st['ratio_par_perp']=1+f['q'];st['frame_beta_if_unshared']=f['df']/winv
Hperp=70/(1+f['q']/2);st['Hperp']=Hperp;st['u_over_c_D46.7']=Hperp*46700/299792.458
json.dump(st,open('fit3.json','w'),indent=1,default=float);print(json.dumps(st,indent=1,default=float))
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':16});fig,ax=plt.subplots(figsize=(12,7))
ax.errorbar(L,Y,yerr=[MLO,MHI],fmt='o',color='k',ms=4,label='Migkas A/A_all')
ll=L.copy()
for name,ff,col in[("full",None,'C0'),("frame0",0.0,'C3')]:
    s=st[name];m=model([s['q'],s['df'],s['l0']],fixf=ff)
    ax.plot(L,m,col,lw=3,label=f"{'best fit' if ff is None else 'best fit, frame speed 0'}: χ²={s['chi2']:.1f}, R²={s['R2']:.2f}")
ax.set_xlabel('Galactic longitude (deg)');ax.set_ylabel('A / A_all')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=1,frameon=False);fig.tight_layout();fig.savefig('hfunc3_fit_migkas.png',dpi=120)
