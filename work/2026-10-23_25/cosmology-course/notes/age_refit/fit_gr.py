"""Exact-GR (LTB + null geodesics) Migkas fit and CMB. usage: fit_gr.py <law> <fwd|back> <gtc|le>"""
import os,sys,json,numpy as np,warnings;warnings.filterwarnings('ignore')
law,mode,vm=sys.argv[1:4]; os.environ['R0DENSE']='1'
if vm=='gtc': os.environ['VMODE']='gtc'
else: os.environ.pop('VMODE',None)
sys.argv=['x',law,mode]
exec(open('fit_lc.py').read().split("ck=f'fwd")[0])
import grltb as G
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
ELLR=np.arange(0,360,4.)
def evg(z):
    try:
        P,o=build(z)
        if o is None: return None
        if vm=='le' and o['U'][o['io']]>F.C: return None   # v<=c case (NM polish unbounded)
        f=G.fields(o); tr=G.trace(f,ELLR,nstep=320,lbmax_frac=0.999); af=G.migkas(tr)   # Migkas grid = fit grid
        if af is None or not np.all(np.isfinite(af)): return None
        chi,l0=F.chi2_l0(af); return P,o,f,tr,af,chi,l0
    except Exception: return None
def cost(z): r=evg(z); return 1e3 if r is None else r[5]
tag=f'{law}_{mode}_{vm}'; ck=f'gr/ck_{tag}.json'
if os.path.exists(ck) and json.load(open(ck)).get('stage')=='nm': x=np.array(json.load(open(ck))['x'])
else:
    r=differential_evolution(cost,BX,seed=11,popsize=8,maxiter=15,polish=False,init='sobol'); json.dump(dict(x=list(map(float,r.x)),stage='de',fun=r.fun),open(ck,'w'))
    r2=minimize(cost,r.x,method='Nelder-Mead',options=dict(maxiter=300,adaptive=True)); x=r2.x if r2.fun<r.fun else r.x
    json.dump(dict(x=list(map(float,x)),stage='nm',fun=float(min(r.fun,r2.fun))),open(ck,'w'))
res=evg(x)
if res is None: json.dump(dict(tag=tag,failed=True),open(f'gr/{tag}.json','w')); sys.exit()
P,o,f,tr,af,chi,l0=res; L=F.local(o); lam=L['Hbar']/H70; io=o['io']
tr=G.trace(f,np.arange(0,360,2.),nstep=500)
A=F.shifted(af,l0)[F.IMIG]; R2=1-np.sum((F.Y-A)**2)/np.sum((F.Y-F.Y.mean())**2)
def nodip(S):
    th=np.radians(np.arange(360.)); X=np.c_[np.ones(360),np.cos(th),np.sin(th)]; c=np.linalg.lstsq(X,np.log(S),rcond=None)[0]; return S/np.exp(c[1]*np.cos(th)+c[2]*np.sin(th))
def orient(S): return np.interp((np.arange(360.)-l0)%360,np.arange(361.),np.append(S,S[0]))
cm={}; curves={}
Su,trec=G.cmb_uniform(f,tr,o); St,tls=G.cmb_Ttraced(f,tr)
for nm_,S in (('uniform_1100',Su),('Ttraced_3000K',St)):
    d=dict(n_finite=int(np.isfinite(S).sum()))
    if np.all(np.isfinite(S)) and np.all(S>0):
        for t_,SS in (('',S),('_nodip',nodip(S))):
            So=orient(SS); r=CC.resid(CC.G_FULL,So); v=[score(SS,p)[0] for p in np.arange(0,360,1.)]
            d['cancel_l0'+t_]=100*(1-CC.rms(r)/S0); d['res_l0'+t_]=CC.rms(r)/1e-5; d['cancel_free'+t_]=100*(1-min(v)/S0)
            d['dTT_modes'+t_]={m['m']:m['A'] for m in CC.fourier(So/So.mean()-1,6)['modes']}; curves[nm_+t_]=(So/So.mean()-1,r)
        d['mean_z']=float(np.mean(1/S-1))
    cm[nm_]=d
rs_D=float(2*o['GM'][io]/F.C**2/o['r'][io])
out=dict(tag=tag,free=nm+['l0'],nfree=len(nm)+1,x=list(map(float,x)),chi2=chi,R2=float(R2),sigma=sig(chi,len(nm)+1),l0=l0,
  M_kg=P['M']*lam,D_Gpc=L['D']*lam,Uo_c=L['Uo']/F.C,rs_over_D=rs_D,Hpar=L['Hpar']*F.KMSMPC/lam,Hper=L['Hper']*F.KMSMPC/lam,age_Gyr=o['t']*lam,
  gk=L['grav_over_kin'],shell_cross=f['cross'],trec_Gyr=trec*lam,cmb=cm,rays_alive_end=int(tr['alive'][-1].sum()))
json.dump(out,open(f'gr/{tag}.json','w'),indent=1,default=float)
fig,ax=plt.subplots(2,1,figsize=(16,10)); ell=np.arange(360.)
ax[0].errorbar(F.IMIG,F.Y,yerr=0.5*(F.MLO+F.MHI),fmt='o',ms=6,c='k',label='Migkas 71 bins');ax[0].plot(ell,F.shifted(af,l0),lw=3,c='C3',label=f'GR model χ²={chi:.2f} σ={out["sigma"]:.1f} l0={l0:.0f}°')
ax[0].set_ylabel('A/A_all',fontsize=18);ax[0].legend(fontsize=14)
for (k,(dtt,r)),c in zip(curves.items(),['C0','C1','C2','C4']):
    ax[1].plot(ell,dtt,lw=3,c=c,label=k+' dT/T'); ax[1].plot(ell,r,lw=2,ls=':',c=c,label=k+f" G·S resid ({cm[k.replace('_nodip','')]['cancel_l0'+('_nodip' if 'nodip' in k else '')]:.1f}%)")
ax[1].plot(ell,CC.G_FULL-1,c='k',lw=1.5,label='G−1');ax[1].set_xlabel('galactic longitude ℓ (deg)',fontsize=18);ax[1].set_ylabel('dT/T, residual',fontsize=18);ax[1].legend(fontsize=10,ncol=2)
if not curves: ax[1].text(0.3,0.5,'CMB: incomplete sky '+str({k:v['n_finite'] for k,v in cm.items()}),transform=ax[1].transAxes,fontsize=16)
for a in ax: a.set_xlim(0,360);a.tick_params(labelsize=14);a.grid()
fig.suptitle(f"EXACT GR {tag}: M={out['M_kg']:.2g} kg D={out['D_Gpc']:.3g} Gpc U_o/c={out['Uo_c']:.3g} r_s/D={rs_D:.2g} H {out['Hpar']:.0f}/{out['Hper']:.0f} age {out['age_Gyr']:.1f}",fontsize=15)
fig.tight_layout();fig.savefig(f'gr/curves_{tag}.png',dpi=80);print('DONE',tag,chi,out['sigma'],flush=True)
