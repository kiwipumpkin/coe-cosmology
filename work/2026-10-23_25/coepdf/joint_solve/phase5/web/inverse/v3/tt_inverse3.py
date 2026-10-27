# FROZEN-FIELD TT for the v3 inverse web (NOT a full refit), exactly as ../tt_inverse.py: model TT = Df (baseline web_RA0.01_M1e22 p8
# best-fit field params & CORR) + Dfix(l=2..128 = Dl of the bg+web residual sky map). Usage: python3 tt_inverse3.py [--dir test/] <seed>...
# Adds D_l at l=20,50,100 vs Planck 2018 TT. Fails on NaN.
import os,sys,json,numpy as np
P5='/workspace/coepdf/joint_solve/phase5/';T='web_RA0.01_M1e22';W=P5+'web/inverse/v3/';av=sys.argv[1:]
if av and av[0]=='--dir':W+=av[1];av=av[2:]
os.environ.update(TAG=T,DG='30.424559192759354',R='21.7439',DFIXF=P5+f'p5_{T}.npz');sys.path.insert(0,P5)
import p8_scratch5 as S
import acoustic_model as M,acoustic_fit as A
from scipy.signal import find_peaks
def fin(x,n):
    if not np.all(np.isfinite(x)):raise SystemExit(f'NONFINITE in {n}')
    return x
ck=json.load(open(P5+f'p8_{T}_ckpt.json'));S.CORR[:]=ck['corr'];b=np.array(json.load(open(P5+f'p8_{T}.json'))['best_x'])
Df=fin(M.Dfield(S.P(b))*S.CORR,'Df');rg=dict(l2_2500=(2,2500),l2_30=(2,30),l30_2500=(30,2500));res={};spec=dict(Lt=A.Lt,Pf=A.Pf,Df=Df)
PL={int(r[0]):float(r[1]) for r in A.Pf};res['planck_D_l20_50_100']={k:PL[k] for k in (20,50,100)}
for sd in av:
  z=np.load(W+f'inv3_s{sd}.npz')
  for nm,key in (('end','Dl_resid'),('start','Dl_start')):
    Dx=np.zeros(M.LMAX+1);Dx[2:129]=z[key][2:129];D=fin(Df+Dx,nm);s={q:A.stats(D,Df,*v) for q,v in rg.items()}
    pk,_=find_peaks(D[100:],prominence=50);pk=pk+100;i=int(np.argmax(D[10:41]))
    res[f's{sd}_{nm}']=dict(chi2_total=s['l2_2500']['chi2_planck'],chi2_2_30=s['l2_30']['chi2_planck'],chi2_30_2500=s['l30_2500']['chi2_planck'],
      r2_total=s['l2_2500']['r2'],r2_2_30=s['l2_30']['r2'],r2_30_2500=s['l30_2500']['r2'],max_l10_40=int(10+i),max_D_l10_40=float(D[10+i]),mean_D_l15_35=float(D[15:36].mean()),
      planck_mean_D_l15_35=float(np.mean([PL[k] for k in range(15,36)])),D_l20=float(D[20]),D_l50=float(D[50]),D_l100=float(D[100]),
      C2=float(D[2]),peaks=[int(p) for p in pk[:5]],peak1_height=float(D[pk[0]]) if len(pk) else None)
    spec[f's{sd}_{nm}']=D;print(sd,nm,res[f's{sd}_{nm}'],flush=True)
json.dump(dict(mode='FROZEN field (baseline p8 params, not refitted)',results=res),open(W+'tt_inverse3.json','w'),indent=1);np.savez(W+'tt_inverse3.npz',**spec)
