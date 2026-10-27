# Phase5 (NEW): table per CMB case (Migkas-profile backgrounds): Migkas chi2/r2, d, k, eta, L, D, mc; bg frac, age, R_ls; TT chi2 total and l=2-30
# (ensemble mean D=Df+Dfix and seed-42 compensated sky), peaks (find_peaks on ensemble D_l, l>=100, prominence>50 uK^2), first-peak height,
# compensation max/rms (p5_fld cancel_err), C2 (D_2 ensemble and seed-42 compensated). Joint = Migkas chi2 + TT chi2 (ensemble, l=2-2500).
import json,glob,os,math,numpy as np
from scipy.signal import find_peaks
O='/workspace/coepdf/joint_solve/phase5/';P=json.load(open('/workspace/coepdf/joint_solve/phase4b/profile_k4b.json'))
KA={round(r[1],2):r[3] for r in P if r[0]=='A'}
def src(t):
    if t.startswith('A_d'):return json.load(open(O+f'fitd_{t[3:]}.json'))
    if t.startswith('K'):return KA[float(t[1:])]
    _,d,lk=t.split('_');return json.load(open(O+f'fitdk_{d}_{lk}.json'))
rows=[]
for line in open(O+'queue6.txt'):
    t=line.split()[0];r=src(t);row=dict(tag=t,d=r['d'],k=r['k'],eta=r['eta'],L=r['L'],D=r['D'],mc=r['mc'],l0=r['l0'],mig_chi2=r['chi2'],r2=r['r2'],mig_age=r['age_Gyr'])
    if os.path.exists(O+f'bg_{t}.json'):
        b=json.load(open(O+f'bg_{t}.json'));row.update({k:b.get(k) for k in ('frac','frac_over_mc','age_Gyr','R_ls','sm_rms','surf','fail')})
    if os.path.exists(O+f'p6_{t}.json'):
        p=json.load(open(O+f'p6_{t}.json'));z=np.load(O+f'p6_{t}.npz');c=json.load(open(O+f'p5_fld_{t}.json'));D=z['D']
        pk,_=find_peaks(D[100:],prominence=50);pk=pk+100
        row.update(TT=p['mean']['l2_2500']['chi2_planck'],TT_2_30=p['mean']['l2_30']['chi2_planck'],TT_comp=p['seed42_comp']['l2_2500']['chi2_planck'],
          TT_comp_2_30=p['seed42_comp']['l2_30']['chi2_planck'],peaks=pk.tolist(),peak1_height=float(D[pk[0]]) if len(pk) else None,
          comp_max=c['cancel_err_max'],comp_rms=c['cancel_err_rms'],D2_ens=p['D2']['ens'],D2_comp=p['D2']['seed42_comp'],comp_frac_at_bound=c['frac_at_bound'],
          joint=r['chi2']+p['mean']['l2_2500']['chi2_planck'])
    rows.append(row)
json.dump(rows,open(O+'summary6.json','w'),indent=1,default=float)
done=[r for r in rows if 'TT' in r]
if done:b=min(done,key=lambda r:r['joint']);json.dump(b,open(O+'best6.json','w'),indent=1,default=float);print('BEST(joint)',b['tag'])
for r in rows:print({k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()})
