# Phase5 web test: TT summary (same definitions as summary6.py) for web_RA0.01_M1e22 -> web/tt_web.json
import json,numpy as np;from scipy.signal import find_peaks
O='/workspace/coepdf/joint_solve/phase5/';t='web_RA0.01_M1e22';p=json.load(open(O+f'p6_{t}.json'));z=np.load(O+f'p6_{t}.npz');c=json.load(open(O+f'p5_fld_{t}.json'));D=z['D']
pk,_=find_peaks(D[100:],prominence=50);pk=pk+100
r=dict(tag=t,TT=p['mean']['l2_2500']['chi2_planck'],TT_2_30=p['mean']['l2_30']['chi2_planck'],TT_comp=p['seed42_comp']['l2_2500']['chi2_planck'],TT_comp_2_30=p['seed42_comp']['l2_30']['chi2_planck'],
 peaks=pk.tolist(),peak1_height=float(D[pk[0]]) if len(pk) else None,D2_ens=p['D2']['ens'],D2_comp=p['D2']['seed42_comp'],comp_max=c['cancel_err_max'],comp_rms=c['cancel_err_rms'],comp_frac_at_bound=c['frac_at_bound'])
json.dump(r,open(O+'web/tt_web.json','w'),indent=1,default=float);print(r)
