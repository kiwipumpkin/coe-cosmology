# v2 diagnostics: (A) t_rec(q) per shell; (B) emission locus r_CoE(t_rec) vs the band of CoE-distances swept by our past light cone; (C) min mismatch per direction.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
z=np.load('cmb_lightpath_v2.npz');D,t0,c,ro=float(z['D']),float(z['t0']),float(z['c']),float(z['ro'])
fig,ax=plt.subplots(1,3,figsize=(16,4.6))
ax[0].semilogx(z['q']*D,z['trec']*t0*1e3);ax[0].axvline(D,ls=':',c='k');ax[0].set_xlabel('shell label q·D [Gpc] (observer shell dotted)');ax[0].set_ylabel('t_rec [Myr]');ax[0].set_title('A. per-shell recombination time (T=3000 K)')
t=np.geomspace(1e-3,1,400);ax[1].fill_between(t*t0,(ro-c*(1-t))*D,(ro+c*(1-t))*D,color='C1',alpha=.4,label='CoE distance of our past light cone')
ax[1].plot(z['trec']*t0,z['xrec']*D,'C0',lw=3,label='emission events (r, t_rec) of all shells')
ax[1].plot(z['tt']*t0,z['robs_t']*D,'k--',label='observer shell r(t)')
ax[1].set_xscale('log');ax[1].set_yscale('log');ax[1].set_xlabel('t [Gyr]');ax[1].set_ylabel('distance from CoE [Gpc]');ax[1].legend(fontsize=7);ax[1].set_title('B. emission locus vs past light cone')
ax[2].plot(np.degrees(z['psi']),z['mismatch']*D);ax[2].set_xlabel('sky angle from anti-CoE [deg]');ax[2].set_ylabel('min |cone - emission event| [Gpc]');ax[2].set_title('C. closest approach per direction')
plt.tight_layout();plt.savefig('cmb_lightpath_v2_trec_locus.png',dpi=110)
