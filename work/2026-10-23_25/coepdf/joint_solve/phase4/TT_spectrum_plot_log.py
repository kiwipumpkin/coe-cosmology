# Phase4 TT power spectrum: CoE model (H0=68.64, age 13.787 Gyr, field refit from scratch) vs Planck 2018
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
z=np.load('p6_age_fld.npz'); l=np.arange(len(z['D']))
b=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-binned_R3.01.txt')
f=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt')
fig,ax=plt.subplots(2,1,figsize=(11,8),sharex=True,gridspec_kw={'height_ratios':[3,1]})
ax[0].plot(f[:,0],f[:,1],'.',ms=1,color='0.75',label='Planck 2018 unbinned')
ax[0].errorbar(b[:,0],b[:,1],yerr=[b[:,2],b[:,3]],fmt='o',ms=3,color='k',label='Planck 2018 binned')
ax[0].plot(l[2:],z['Lt'][2:],'b-',lw=1.2,label=r'$\Lambda$CDM best fit ($\chi^2$=2558)')
ax[0].plot(l[2:],z['D'][2:],'r-',lw=1.5,label=r'CoE ensemble mean ($\chi^2\approx$12,471)')
ax[0].plot(l[2:],z['Drc'][2:],color='orange',lw=0.6,alpha=0.7,label=r'CoE seed 42, compensated ($\chi^2\approx$13,029)')
ax[0].set_ylabel(r'$D_\ell$ [$\mu$K$^2$]'); ax[0].set_ylim(0,7000); ax[0].legend(fontsize=9)
ax[0].set_title(r'TT power spectrum, $\langle H_0\rangle$=68.64, age 13.787 Gyr, $R_{ls}$=21.77 Gpc (phase 4), log $\ell$')
for p in [224,514,822,1124,1405]: ax[0].axvline(p,color='k',ls=':',lw=0.5)
lb=np.r_[f[:28,0],b[:,0]].astype(int); bd=np.r_[f[:28,1],b[:,1]]; err=np.r_[0.5*(f[:28,2]+f[:28,3]),0.5*(b[:,2]+b[:,3])]
ax[1].plot(lb,(np.interp(lb,l,z['D'])-bd)/err,'r.',label='CoE')
ax[1].plot(lb,(np.interp(lb,l,z['Lt'])-bd)/err,'b.',label=r'$\Lambda$CDM')
ax[1].axhline(0,color='k',lw=0.5); ax[1].set_ylabel(r'(model$-$data)/$\sigma$'); ax[1].set_xlabel(r'$\ell$'); ax[1].legend(fontsize=8)
ax[1].set_xlim(2,2500); ax[1].set_xscale('log'); plt.tight_layout(); plt.savefig('TT_spectrum_plot_log.png',dpi=140)
