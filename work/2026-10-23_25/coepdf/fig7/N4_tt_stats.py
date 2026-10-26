# NEW (rev7): fit statistics and peak table for the current TT solve (p6_p8) vs LCDM, one peak finder (Savitzky-Golay 61/3, prominence 20, min separation 200).
import os,sys,json,numpy as np
from scipy.signal import savgol_filter,find_peaks
P='/workspace/coepdf/joint_solve/phase2/';z=np.load(P+'p6_p8.npz');J=json.load(open(P+'p6_p8.json'))
os.chdir('/workspace/redweb2/stochfield');sys.path.insert(0,'.');import acoustic_fit as A
Lt=z['Lt'];out={'lcdm':{k:A.stats(Lt,0*Lt,*r) for k,r in dict(l2_30=(2,30),l30_2500=(30,2500),l2_2500=(2,2500)).items()}}
out['lcdm_fit_objective']=float(A.chi2_fit(Lt));out['model_fit_objective']=float(A.chi2_fit(z['D']))
N=2499;pm=15;pl=6
for k in ('mean','seed42','seed42_comp'):
    c=J[k]['l2_2500']['chi2_planck'];c0=out['lcdm']['l2_2500']['chi2_planck']
    out[k]=dict(chi2=c,dchi2=c-c0,dBIC=(c-c0)+(pm-pl)*np.log(N),dAIC=(c-c0)+2*(pm-pl))
Pf=z['Pf'];l=Pf[:,0].astype(int)
def peaks(y,ls):
    s=savgol_filter(y,61,3);i,_=find_peaks(s,prominence=20,distance=200);i=i[(ls[i]>100)][:5];return [(int(ls[j]),float(round(s[j]))) for j in i]
L=np.arange(len(z['D']))
out['peaks']=dict(planck=peaks(Pf[:,1],l),model_mean=peaks(z['D'][2:2501],L[2:2501]),lcdm=peaks(Lt[2:2501],L[2:2501]))
json.dump(out,open(__file__.replace('.py','.json'),'w'),indent=1,default=float);print(json.dumps(out,default=float))
