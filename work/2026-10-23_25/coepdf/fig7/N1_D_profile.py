# NEW in rev7: choice of background. Migkas chi2 (profile refit, eta/mc/l0 free, CatWISE lock) vs D (p1) and the
# LTB-implied H0, age and R_ls vs D (p3). Bins are correlated: chi2 differences are descriptive.
import json,glob,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
P='/workspace/coepdf/joint_solve/phase2/';r=[list(json.load(open(f)).values())[0] for f in glob.glob(P+'p1_migkas_D*.json')];r.sort(key=lambda v:v['D'])
B=json.load(open(P+'p3_background_relations.json'))['by_D'];D=[v['D'] for v in r]
fig,ax=plt.subplots(1,2,figsize=(11,4))
ax[0].plot(D,[v['chi2'] for v in r],'ko-');ax[0].axhline(9.27823,color='r',ls='--',label='locked fit, $\\chi^2$=9.278');ax[0].axvline(34.96,color='m',ls=':',label='chosen D=34.96 Gpc')
ax[0].set_xlabel('D [Gpc]');ax[0].set_ylabel('Migkas $\\chi^2$ (71 correlated cones)');ax[0].legend(fontsize=8)
k=sorted(B,key=float);Dk=[float(x) for x in k];a2=ax[1].twinx()
ax[1].plot(Dk,[B[x]['H0_ltb_kms_Mpc'] for x in k],'bo-',label='$H_0$ (LTB) [km/s/Mpc]');a2.plot(Dk,[B[x]['age_ltb_Gyr'] for x in k],'rs-',label='age (LTB) [Gyr]')
ax[1].set_xlabel('D [Gpc]');ax[1].set_ylabel('$H_0$ [km/s/Mpc]',color='b');a2.set_ylabel('age [Gyr]',color='r');ax[1].axvline(34.96,color='m',ls=':')
fig.tight_layout();fig.savefig('N1_D_profile.png',dpi=140,facecolor='white')
