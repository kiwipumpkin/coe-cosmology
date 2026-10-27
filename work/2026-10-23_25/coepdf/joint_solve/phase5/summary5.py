# Phase5 (NEW): table per d: Migkas chi2/r2 (fitd_*.json), TT chi2 total & l=2-30 (p6, ensemble mean and seed-42 compensated sky),
# peaks (scipy find_peaks on ensemble D_l, l>=100, prominence>50 uK^2), first-peak height, compensation max/rms (p5_fld cancel_err), D_2.
import json,glob,os,numpy as np
from scipy.signal import find_peaks
O='/workspace/coepdf/joint_solve/phase5/';rows=[]
for f in sorted(glob.glob(O+'fitd_*[0-9].json'),key=lambda f:json.load(open(f))['d']):
    r=json.load(open(f));T='D%.2f'%(round(r['D']*20)/20);row=dict(d=r['d'],D=r['D'],mc=r['mc'],l0=r['l0'],eta=r['eta'],k=r['k'],L=r['L'],age=r['age_Gyr'],mig_chi2=r['chi2'],r2=r['r2'],flags=r['flags'],tag=T)
    if os.path.exists(O+f'p6_{T}.json'):
        p=json.load(open(O+f'p6_{T}.json'));z=np.load(O+f'p6_{T}.npz');c=json.load(open(O+f'p5_fld_{T}.json'));lt=json.load(open(O+f'ltb_scan_{T}.json'))[0]
        D=z['D'];pk,_=find_peaks(D[100:],prominence=50);pk=pk+100
        row.update(R_ls=lt['eta_rec_D']*r['D'],TT_mean=p['mean']['l2_2500']['chi2_planck'],TT_mean_2_30=p['mean']['l2_30']['chi2_planck'],
          TT_comp=p['seed42_comp']['l2_2500']['chi2_planck'],TT_comp_2_30=p['seed42_comp']['l2_30']['chi2_planck'],peaks=pk.tolist(),
          peak1_height=float(D[pk[0]]) if len(pk) else None,comp_max=c['cancel_err_max'],comp_rms=c['cancel_err_rms'],D2_ens=p['D2']['ens'],D2_comp=p['D2']['seed42_comp'],
          comp_frac_at_bound=c['frac_at_bound'])
    rows.append(row)
json.dump(rows,open(O+'summary5.json','w'),indent=1)
done=[r for r in rows if 'TT_mean' in r]
if done:b=min(done,key=lambda r:r['TT_mean']);json.dump(b,open(O+'best_d.json','w'),indent=1);print('BEST',b)
for r in rows:print({k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()})
