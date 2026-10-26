# rev6 Fig 20 (low multipoles) regenerated from the rev7 solve: p6_p8 (D=34.96 Gpc LTB background, p8 field, p5_c_p8 bodies).
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase2/';Z=np.load(O+'p6_p8.npz');J=json.load(open(O+'p6_p8.json'))
Pf,Lt=Z['Pf'],Z['Lt'];ll=np.arange(len(Z['D']));r=slice(2,31);s=Pf[:,0]<=30
fig,a=plt.subplots(figsize=(7,3.6))
a.errorbar(Pf[s,0],Pf[s,1],[Pf[s,2],Pf[s,3]],fmt='o',color='grey',ms=3,label='Planck')
a.plot(ll[r],Lt[r],'k',label='ΛCDM');a.plot(ll[r],Z['D'][r],'r',label='CoE ensemble mean')
a.plot(ll[r],Z['Drc'][r],'m--',label='CoE one-sky (seed 42 + bodies)')
a.set_xlabel('ℓ');a.set_ylabel('D_ℓ (μK²)');a.set_xlim(1.5,30.5);a.legend(fontsize=7);fig.tight_layout()
fig.savefig(__file__.replace('.py','.png'),dpi=130);print('ok')
