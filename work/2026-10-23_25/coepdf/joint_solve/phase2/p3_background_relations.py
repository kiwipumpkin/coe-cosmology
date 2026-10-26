# Phase2 A3: implied H0, age, last-scattering comoving distance for the chosen D, from the SAME LTB used to trace the CMB background
# (famP/ltb.py poly C=[-0.041,-0.035,0.1353,0], frac 20000, u=4; units D=1, c=1 so time unit D/c, H unit c/D),
# and from the Migkas prof_core model (age, H t reported by p1).
import json,numpy as np
C=299792.458;GYR_PER_GPC=3.26156
T=np.load('/workspace/redweb2/realweb/los_flow_R/tab.npz');L=np.load('/workspace/redweb2/realweb/los_flow/ltb.npz')
t0=float(T['t0']);eta_rec=float(T['eta'][0]);H_loc=float(L['H'][-1]);Om=float(L['Omega'][-1])
out={'ltb':dict(t0_Dc=t0,eta_rec_D=eta_rec,H_local_cD=H_loc,H_t0=H_loc*t0,Omega_local_today=Om),'by_D':{}}
p1={}
for D in (30,34.96,42,50,58.2):
    try:p1[D]=json.load(open(f'/workspace/coepdf/joint_solve/phase2/p1_migkas_D{D:g}.json'))[str(float(D))]
    except Exception:pass
for D in (34.96,42,58.2):
    o=dict(D=D,H0_ltb_kms_Mpc=H_loc*C/(D*1000),age_ltb_Gyr=t0*D*GYR_PER_GPC,R_ls_ltb_Gpc=eta_rec*D,
           H0_needed_for_13p8Gyr_given_Ht=H_loc*t0*977.8/13.8)
    if D in p1:o.update(migkas_chi2=p1[D]['chi2'],migkas_r2=p1[D]['r2_pearson'],migkas_age_Gyr_profcore=p1[D]['age_Gyr'],migkas_Ht_profcore=p1[D]['Ht'],
                         H0_profcore_if_13p8Gyr=p1[D]['Ht']*977.8/13.8)
    out['by_D'][str(D)]=o
json.dump(out,open('/workspace/coepdf/joint_solve/phase2/p3_background_relations.json','w'),indent=1);print(json.dumps(out,indent=1))
