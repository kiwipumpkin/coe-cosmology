# v3C: quadrupole/octopole axes (max angular-momentum dispersion) of the monopole+dipole removed map vs AoE (260,60); writes cmb_lightpath_v3C.json (all key numbers).
import numpy as np,healpy as hp,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from cmb_lightpath_v3_common import *
m,m1,cl,cl1,kind,mean=stats('C');ax2,ax3=axis(m1,2),axis(m1,3);l,pm,pp=pct(cl1);ll=np.arange(len(cl1))
J=dict(rule='C',kind=kind,mean_profile_value=mean,rms_uK_monopole_removed=float(np.std(m)*T0uK),rms_uK_dipole_removed=float(np.std(m1)*T0uK),
 dipole_amp_dTT=float(np.sqrt(3*cl[1]/(4*np.pi))),Dl_uK2_l2_3_4_10=[float(ll[i]*(ll[i]+1)/2/np.pi*cl1[i]*T0uK**2) for i in (2,3,4,10)],
 pct_l2_30_model=pm.tolist(),pct_l2_30_planck=pp.tolist(),axis_l2=ax2,axis_l3=ax3)
json.dump(J,open('cmb_lightpath_v3C.json','w'),indent=1);print(json.dumps({k:v for k,v in J.items() if not k.startswith('pct')}))
fig=plt.figure(figsize=(9,5));hp.mollview(m1*T0uK,fig=fig.number,title='v3C monopole+dipole removed [uK], l=2 (green), l=3 (magenta) axes',cmap='coolwarm')
for a,cc in((ax2,'g'),(ax3,'m')):hp.projscatter([a['l'],a['l']+180],[a['b'],-a['b']],lonlat=True,c=cc,s=80)
hp.projscatter(260,60,lonlat=True,c='k',marker='x',s=100);hp.projscatter(180,-15,lonlat=True,c='purple',marker='*',s=120)
plt.savefig('cmb_lightpath_v3C_lowl_axes.png',dpi=110,bbox_inches='tight')
