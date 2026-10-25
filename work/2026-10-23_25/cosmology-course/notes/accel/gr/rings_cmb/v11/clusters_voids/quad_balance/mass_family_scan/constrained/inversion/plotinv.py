import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from obs import *;R=json.load(open('inv1_results.json'));Ga,_=obs_grid('/workspace/inv/smica_alm_l30.fits');L=np.arange(2,31)
for n in['red_inv_t0.5','refpoly_inv_t1.0']:
    r=R[n];i=hp.ang2pix(NA,*r['axis_Agfree'],lonlat=True);am=np.array(r['a_motion']);ag=np.array(r['a_grav'])*r['Ag_free']
    f,ax=plt.subplots(figsize=(9,5.5));ax.plot(L,am[2:],'o-',label='motion term');ax.plot(L,ag[2:],'s-',label='gravity term ×A_g=%.3g'%r['Ag_free']);ax.plot(L,(am+ag)[2:],'k--',label='sum')
    ax.plot(L,Ga[2:,i],'r^-',label='Planck a_l0, axis (%.1f,%.1f), chi2=%.2f'%(*r['axis_Agfree'],r['chi2_Agfree']));ax.fill_between(L,-np.sqrt(ClP[2:]),np.sqrt(ClP[2:]),color='gray',alpha=.2)
    ax.set_yscale('symlog',linthresh=1e-1);ax.set_xlabel('ℓ');ax.set_ylabel('a_ℓ0 [μK]');ax.set_title(n);ax.legend(fontsize=7);f.tight_layout();f.savefig(f'al0_{n}.png',dpi=120)
    print(n,np.round(r['Dl_motion'][2:7],1),max(r['Dl_motion'][7:]))
