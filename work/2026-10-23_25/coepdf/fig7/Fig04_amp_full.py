# REGENERATED for rev7: full Migkas A/A_all profile with the locked CoE eta-edge curve (priority-1 fit, unchanged),
# the D=34.96 Gpc profile refit (p1) and the locked curve after the phase-2 compensating bodies (p7 method, final run p5_c_p8).
# chi2 uses the mean of the asymmetric errors over the 71 overlapping cones (descriptive only; bins are correlated).
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high'];s=0.5*(el+eh)
c=np.loadtxt('/workspace/cosmology-course/notes/catwise-eta-migkas-fit-d14.3.csv',delimiter=',',skiprows=1)
P1=json.load(open('/workspace/coepdf/joint_solve/phase2/p1_migkas_D34.96.json'))['34.96'];cf=np.array(P1['curve']);lf=np.arange(len(cf))*360/len(cf)
B=json.load(open('/workspace/coepdf/joint_solve/phase2/p7_migkas_bodies.json'))['p5_c_p8.npz']
m0=np.interp(l,c[:,0],c[:,1],period=360);m1=np.interp(l,lf,cf,period=360);chi=lambda m:np.sum((y-m)**2/s**2)
fig,ax=plt.subplots(2,1,figsize=(12,6.4),sharex=True,gridspec_kw=dict(height_ratios=[3,1.2]))
ax[0].errorbar(l,y,[el,eh],fmt='ko',ms=3.5,capsize=2,label='Migkas et al. (2020) A/A$_{\\rm all}$, 71 sliding 90° cones')
ax[0].plot(c[:,0],c[:,1],'r',lw=2.2,label=f'locked CoE $\\eta$-edge fit (priority 1), $\\chi^2$={chi(m0):.3f}')
ax[0].plot(lf,cf,'b--',lw=1.4,label=f'profile refit at D=34.96 Gpc, $\\chi^2$={P1["chi2"]:.3f}')
ax[0].plot([],[],' ',label=f'with compensating bodies: $\\chi^2$={B["chi2_after"]:.5f}, max|$\\Delta A$|={B["max_dA"]:.1e} (curve indistinguishable)')
ax[0].set_ylabel('A/A$_{\\rm all}$');ax[0].legend(fontsize=8.5,loc='upper right');ax[0].grid(alpha=.3)
ax[1].errorbar(l,y-m0,[el,eh],fmt='ko',ms=3,capsize=2);ax[1].axhline(0,color='r');ax[1].set_ylabel('data − locked');ax[1].set_xlabel('Galactic longitude [deg]');ax[1].set_xlim(0,360);ax[1].grid(alpha=.3)
fig.tight_layout();fig.savefig('Fig04_amp_full.png',dpi=140,facecolor='white');print(chi(m0),chi(m1))
