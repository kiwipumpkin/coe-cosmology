# Phase5 web TT ablation (frozen stochastic field): same pipeline pieces as tt_web.sh (p5_comp5 kernel/Dl, p8 Dfield best-fit params+CORR, p6 stats),
# baseline tag web_RA0.01_M1e22. Model TT = Df(field, best-fit params of the baseline p8 fit) + Dfix(l=2..128 = Dl of the bg+bodies residual map).
# Variants: base | noweb (Dfix from bg map B alone) | noISW (field recomputed with ISW kernel KI=0; Dfix baseline) |
#           nonear (nearest 10% of bodies by |y| removed from residual) | obsonly (diagnostic: Dfix = Dl of target observed sky Tobs).
# Field parameters are NOT refitted per variant (isolates the cause). Fails loudly on non-finite values.
import os,sys,json,numpy as np,healpy as hp
P5='/workspace/coepdf/joint_solve/phase5/';T='web_RA0.01_M1e22';OUT=P5+'web/ablation/'
os.environ.update(TAG=T,DG='30.424559192759354',R='21.7439',DFIXF=P5+f'p5_{T}.npz');sys.path.insert(0,P5)
import p8_scratch5 as S   # builds M.EX at the baseline R with the p8 Pd; chdir to stochfield
import acoustic_model as M,acoustic_fit as A
from scipy.signal import find_peaks
def fin(x,n):
    if not np.all(np.isfinite(x)):raise SystemExit(f'NONFINITE in {n}')
    return x
ck=json.load(open(P5+f'p8_{T}_ckpt.json'));S.CORR[:]=ck['corr'];b=np.array(json.load(open(P5+f'p8_{T}.json'))['best_x'])
z=np.load(P5+f'p5_{T}.npz');y,Mb,B,Resid,Tobs=z['y'],z['M'],z['B'],z['Resid'],z['Tobs']
NS=64;T0=2.7255e6;R=21.7439;EPS=0.025*R;MSUN=1.989e30;KG=7.7e51/0.0039;RV=14.3;n=np.array(hp.pix2vec(NS,np.arange(hp.nside2npix(NS)))).T
def kern(v):
    d=v[None,:]-R*n;dd=(np.sum(d*d,1)+EPS**2)**1.5;k=np.sum(n*d,1)/dd-n@v/np.linalg.norm(v)**3;return (k-k.mean())*MSUN/KG*RV**2
def Dl(m,L=3*NS-1):l=np.arange(L+1);return l*(l+1)/(2*np.pi)*hp.anafast(hp.remove_dipole(m),lmax=L)*T0**2
# sanity: residual reproduces from kernel
chk=float(np.abs(B+sum(kern(y[i])*Mb[i] for i in range(len(Mb)))-Resid).max());print('kernel reproduce max err',chk)
assert chk<1e-6*max(1,np.abs(Resid).max())+1e-9
r=np.linalg.norm(y,axis=1);near=r<=np.quantile(r,0.10)
ang=np.degrees(2*np.arctan(EPS/r))   # angular size of the softening core seen from observer
Rnear=Resid-sum(kern(y[i])*Mb[i] for i in np.where(near)[0])
def dfix(m):d=np.zeros(M.LMAX+1);a=Dl(m);d[2:129]=a[2:129];return fin(d,'Dfix')
Df0=fin(M.Dfield(S.P(b))*S.CORR,'Df')
KI0=M.EX.KI.copy();M.EX.KI=np.zeros_like(KI0);Df_noisw=fin(M.Dfield(S.P(b))*S.CORR,'Df_noISW');M.EX.KI=KI0
V=dict(base=(Df0,dfix(Resid)),noweb=(Df0,dfix(B)),noISW=(Df_noisw,dfix(Resid)),nonear=(Df0,dfix(Rnear)),obsonly=(Df0,dfix(Tobs)))
rg=dict(l2_2500=(2,2500),l2_30=(2,30),l30_2500=(30,2500));res={};spec={}
for k,(Df,Dx) in V.items():
    D=fin(Df+Dx,k);spec[k]=D;s={q:A.stats(D,Df,*v) for q,v in rg.items()}
    lw=np.arange(10,41);i=int(np.argmax(D[10:41]));pk,_=find_peaks(D[100:],prominence=50);pk=pk+100
    res[k]=dict(chi2_total=s['l2_2500']['chi2_planck'],chi2_2_30=s['l2_30']['chi2_planck'],chi2_30_2500=s['l30_2500']['chi2_planck'],
      r2_total=s['l2_2500']['r2'],r2_2_30=s['l2_30']['r2'],bump_l=int(lw[i]),bump_height=float(D[10+i]),C2=float(D[2]),
      first_peak_l=int(pk[0]),first_peak_height=float(D[pk[0]]),Dfix_20=float(Dx[20]),Df_20=float(Df[20]))
    print(k,{q:(round(v,3) if isinstance(v,float) else v) for q,v in res[k].items()},flush=True)
meta=dict(tag=T,mode='frozen field (baseline p8 best-fit params & CORR, not refitted)',near_def='nearest 10% of bodies by distance from observer',
  n_near=int(near.sum()),r_near_max_Gpc=float(r[near].max()),R_Gpc=R,n_bodies=len(Mb),n_near_with_core_gt5deg=int((ang[near]>5).sum()),n_all_core_gt5deg=int((ang>5).sum()),
  near_mass_pos=float(Mb[near][Mb[near]>0].sum()),near_mass_neg=float(Mb[near][Mb[near]<0].sum()),kernel_reproduce_err=chk,
  ISW_switch='acoustic_model.EX.KI set to 0 in-memory (no file edits)',baseline_tt_web_json=json.load(open(P5+'web/tt_web.json')))
json.dump(dict(meta=meta,variants=res),open(OUT+'ablation_metrics.json','w'),indent=1,default=float)
np.savez(OUT+'ablation_spectra.npz',Lt=A.Lt,Pf=A.Pf,**spec,Dl_resid_base=Dl(Resid),Dl_B=Dl(B),Dl_nonear=Dl(Rnear),Dl_obs=Dl(Tobs))
print('META',meta)
