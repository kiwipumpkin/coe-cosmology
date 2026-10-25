"""rec_fit7: rec_fit6 variant A + TD only; Migkas-only objective; CoE pinned at l=180 (shifted(af,180)); comparison with l0 free."""
import json,numpy as np,sys;sys.path.insert(0,'.')
exec(open('rec_fit6.py').read().split('WT,WC=')[0])
from bubble_core import shifted,IMIG
VAR=('A',1)
def c2fix(af): return float(np.sum(F.W*(shifted(af,180.)[IMIG]-F.Y)**2))
def R2(af,l0):
    m=shifted(af,l0)[IMIG]; yb=np.sum(F.W*F.Y)/np.sum(F.W); return float(1-np.sum(F.W*(m-F.Y)**2)/np.sum(F.W*(F.Y-yb)**2))
def obj(p,mode):
    try:
        af=migkas(p,VAR,6)
        if af is None or np.any(af<=0) or not np.all(np.isfinite(af)): return 1e6
        return c2fix(af) if mode=='pinned' else F.chi2_l0(af)[0]
    except Exception: return 1e6
BOX=[(-3,4),(-6,6),(np.log(10),np.log(3000))]; out={}
for mode in ('pinned','free'):
    r=differential_evolution(obj,BOX,args=(mode,),seed=7,popsize=15,maxiter=80,tol=1e-9,polish=False)
    r2=minimize(obj,r.x,args=(mode,),method='Nelder-Mead',options=dict(maxiter=4000)); x=r2.x if r2.fun<r.fun else r.x
    res=dict(x=x.tolist())
    for st in (6,1):
        ell=np.arange(0,360,float(st)); q=terms(x,ell,None,VAR); af=migkas(x,VAR,st)
        c2,l0=(c2fix(af),180.) if mode=='pinned' else F.chi2_l0(af); T=3000/q['onez']; h=harm(T,ell)
        res[f'{st}deg']=dict(chi2=c2,sigma=sig(c2,6),l0=l0,R2=R2(af,l0),Tmean=float(T.mean()),dipole=h['dip'],quad=h['quad'],
          rms_nodip_frac=h['rms_nodip'],rms_nodip_K=h['rms_nodip']*float(T.mean()),rms_nodip_1e5=h['rms_nodip']/1e-5,
          maxabs_nodip_K=float(np.abs(h['r']).max()*T.mean()),af=af.tolist(),
          terms={t:[float(q[t].min()),float(q[t].mean()),float(q[t].max())] for t in ('stretch','grav','gamma','dop','onez')})
    res['params']={t:float(q[t]) for t in ('U0','mu','ts','D','U','H','R')}; res['params']['2mu_over_min_re']=float(2*q['mu']/q['r_e'].min())
    out[mode]=res; print(mode,res['params'],{s:{k:res[s][k] for k in ('chi2','sigma','l0','R2','Tmean','dipole','quad','rms_nodip_K','rms_nodip_1e5')} for s in ('6deg','1deg')},flush=True)
json.dump(out,open('rec_fit7.json','w'),indent=1,default=float)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':22})
fig,ax=plt.subplots(2,1,figsize=(22,18));ang=np.arange(360.)
for a,mode in zip(ax,('pinned','free')):
    r=out[mode]['1deg'];a.errorbar(IMIG,F.Y,yerr=1/np.sqrt(F.W),fmt='ko',ms=6,label='Migkas A/A_all (71 bins)')
    a.plot(ang,shifted(np.array(r['af']),r['l0']),lw=4,label=f"model, CoE at ℓ={r['l0']:.0f}°");a.axvline(r['l0'],ls='--',c='C3')
    a.set_title(f"rec_fit7 {mode}: χ²={r['chi2']:.1f} (null 45.1), σ={r['sigma']:.2f}, R²={r['R2']:.2f}; ⟨T⟩={r['Tmean']:.3g} K",fontsize=24);a.set_xlabel('ℓ [deg]');a.grid();a.legend()
fig.tight_layout();fig.savefig('rec_fit7_migkas.png',dpi=60)
