"""Corrected light-cone Migkas fits. usage: fit_lc.py <law> <fwd|back>"""
import sys,os,json,numpy as np,warnings;warnings.filterwarnings('ignore')
law,mode=sys.argv[1],sys.argv[2]; GTC=os.environ.get('VMODE')=='gtc'; UFIX=float(os.environ.get('UFIX','0')); MAXIT=int(os.environ.get('MAXIT','25'))
TAG=('_gtc' if GTC else '')+(f'_U{UFIX:g}' if UFIX else ''); os.environ['VMAX']='1'; sys.argv=['x',law]
exec(open('sr_local.py').read().split("if __name__")[0])
from scipy.optimize import differential_evolution,minimize
SH=law in('eta','pl','soft')
if mode=='fwd': nm=list(names); BX=[B[n] for n in nm]
else:
    nm=(['vo'] if not GTC else (['lv'] if not UFIX else []))+['p','lM','lf','lL']+(['s'] if SH else [])
    BXd={'vo':(1e-3,1.0),'lv':(0,12),'p':(-3,4),'lM':(40,60),'lf':(-4,4),'lL':(-3,5),'s':(0.2,6)}; BX=[BXd[n] for n in nm]   # lv=log10(U_o/c): search box only
def build(z):
    if mode=='fwd':
        P=P_of(z); o=F.integrate(P,True)
        Uo=o['U'][o['io']] if o is not None else 0
        if o is None or (Uo>F.C if not GTC else Uo<=F.C): return None,None
        return P,o
    d=dict(zip(nm,z)); M=10.0**d['lM']
    if GTC: d['vo']=UFIX if UFIX else 10**d['lv']
    if GTC and d['vo']<=1: return None,None   # Chris's U_o>c (NM polish is unbounded)
    P=dict(kind=law,M=M,L=10**d['lL'],p=d['p'],s=d.get('s'),rhoa=1.0)
    rho,dm,Me=F.setup(P); io=np.searchsorted(F.R0,1.0); ms=np.sum(dm[:io])+dm[io]/2; P['rhoa']=10**d['lf']*M/ms
    P['Uo']=d['vo']*F.C; P['f_matter']=10**d['lf']; o=F.integrate_back(P); return P,o
def evl(z):
    try:
        P,o=build(z)
        if o is None: return None
        af=F.migkas_lightcone(o)
        if af is None or not np.all(np.isfinite(af)): return None
        chi,l0=F.chi2_l0(af); return P,o,af,chi,l0
    except Exception: return None
def cost(z):
    r=evl(z)
    if r is None: return 1e3
    pen=0.0
    if UFIX and mode=='fwd': pen=100*(np.log10(r[1]['U'][r[1]['io']]/F.C)-np.log10(UFIX))**2   # profile penalty (scan only)
    return r[3]+pen
ck=f'fwd/lc_{law}_{mode}{TAG}_ckpt.json'
if os.path.exists(ck): x=np.array(json.load(open(ck))['x'])
else:
    r=differential_evolution(cost,BX,seed=7,popsize=8,maxiter=MAXIT,polish=False,init='sobol'); json.dump(dict(x=list(map(float,r.x)),stage='de'),open(ck,'w'))
    r2=minimize(cost,r.x,method='Nelder-Mead',options=dict(maxiter=500,adaptive=True)); x=r2.x if r2.fun<r.fun else r.x
    json.dump(dict(x=list(map(float,x)),stage='nm'),open(ck,'w'))
P,o,af,chi,l0=evl(x); L=F.local(o); lam=L['Hbar']/H70; io=o['io']
A=F.shifted(af,l0)[F.IMIG]; R2=1-np.sum((F.Y-A)**2)/np.sum((F.Y-F.Y.mean())**2)
afold,_=F.migkas_forward(o); chi_old,_=F.chi2_l0(afold)
def nodip(S):
    th=np.radians(np.arange(360.)); X=np.c_[np.ones(360),np.cos(th),np.sin(th)]; c=np.linalg.lstsq(X,np.log(S),rcond=None)[0]
    return S/np.exp(c[1]*np.cos(th)+c[2]*np.sin(th))
cm={}
for q in (0.1,0.23,0.5,1.0):
    c=F.cmb_sr_local(o,q); S=c['S']
    if not np.all(np.isfinite(S)): cm[q]=dict(nfound=int(np.isfinite(S).sum())); continue
    r={}
    for tag,SS in (('',S),('_nodip',nodip(S))):
        rl,res=score(SS,l0); v=[score(SS,p)[0] for p in np.arange(0,360,1.)]; pb=float(np.argmin(v))
        r['cancel_l0'+tag]=100*(1-rl/S0); r['res_l0'+tag]=rl/1e-5; r['cancel_free'+tag]=100*(1-min(v)/S0); r['res_free'+tag]=min(v)/1e-5; r['phi'+tag]=pb
        r['modes'+tag]={m['m']:m['A'] for m in CC.fourier(res,6)['modes']}
    r['dtLS_0m180']=float((c['te'][0]-c['te'][180])*lam); cm[q]=r
Ht=[]; 
for frac in (1.0,0.5,0.1,0.01):
    k=int(np.argmin(abs(o['ts']-frac*o['t']))); Ht.append((frac,float(o['Us'][k][io]/o['rs'][k][io]*F.KMSMPC/lam)))
out=dict(law=law,mode=mode,free=nm+['l0'],nfree=len(nm)+1,x=list(map(float,x)),chi2=chi,R2=float(R2),sigma=sig(chi,len(nm)+1),l0=l0,
  chi2_old_todaypositions_same_params=chi_old,M_kg=P['M']*lam,D_Gpc=L['D']*lam,vo_c=L['Uo']/F.C,Hpar=L['Hpar']*F.KMSMPC/lam,Hper=L['Hper']*F.KMSMPC/lam,
  age_Gyr=o['t']*lam,gk=L['grav_over_kin'],rs_over_D=float(2*o['GM'][io]/F.C**2/L['D']),p=P['p'],L_Gpc=P['L']*lam,shape=P['s'],f=P.get('f_matter'),
  H_obs_vs_t=Ht,cmb=cm,shell_cross=o['cross'])
uni={}
for ar in (1/1100.,1/100.):
    c=F.cmb_uniform(o,ar); S=c['S']; d=dict(trec_Gyr=c['trec']*lam,lookback_Gyr=c['lookback']*lam,mean_z=c['mean_z'],n_inside=int(c['inside'].sum()))
    if np.all(np.isfinite(S)):
        for tag,SS in (('',S),('_nodip',nodip(S))):
            rl,res=score(SS,l0); v=[score(SS,p)[0] for p in np.arange(0,360,1.)]
            d['cancel_l0'+tag]=100*(1-rl/S0); d['res_l0'+tag]=rl/1e-5; d['cancel_free'+tag]=100*(1-min(v)/S0); d['phi'+tag]=float(np.argmin(v))
            d['S_modes'+tag]={m['m']:m['A'] for m in CC.fourier(SS/SS.mean()-1,6)['modes']}; d['res_modes'+tag]={m['m']:m['A'] for m in CC.fourier(res,6)['modes']}
        d['T_ell_K']=(2.725*S/S.mean()).tolist()
    uni[str(round(1/ar))]=d
out['cmb_uniform']=uni
if mode=='back':
    tb=o['tb']*lam; sel=(F.R0>0.1)&(F.R0<10)
    out.update(reached_origin=o['reached_origin'],age_kepler=o['age_kepler']*lam,bang_time_spread_Gyr=float(np.nanmax(tb[sel])-np.nanmin(tb[sel])),bang_time_range=[float(np.nanmin(tb[sel])),float(np.nanmax(tb[sel]))])
json.dump(out,open(f'fwd/lc_{law}_{mode}{TAG}'+('_gr' if os.environ.get('GRZ')=='gr' else '')+'.json','w'),indent=1,default=float)
print('DONE',{k:v for k,v in out.items() if k!='cmb'},flush=True)
