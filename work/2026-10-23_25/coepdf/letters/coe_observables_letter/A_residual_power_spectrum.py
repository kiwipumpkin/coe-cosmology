# Power spectrum P_m = amp^2 of the A/A_all residual harmonics, from A_fit_residual_fourier.json
import json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=json.load(open('A_fit_residual_fourier.json'));M=d['modes']
m=[x['m'] for x in M];P=[x['amp']**2 for x in M];f=[x['frac_power'] for x in M]
fig,ax=plt.subplots(1,2,figsize=(11,4))
b=ax[0].bar(m,P,color='steelblue')
for i,(pp,ff) in enumerate(zip(P,f)):ax[0].text(m[i],pp,'%.1f%%'%(100*ff),ha='center',va='bottom',fontsize=8)
ax[0].set(xlabel='harmonic m (in Galactic longitude)',ylabel='power  amp$^2$',title='Residual power spectrum (linear)')
ax[1].semilogy(m,P,'o-',color='steelblue');ax[1].set(xlabel='harmonic m',ylabel='power amp$^2$ (log)',title='Residual power spectrum (log)')
for a in ax:a.set_xticks(m);a.grid(alpha=.3)
plt.tight_layout();plt.savefig('A_residual_power_spectrum.png',dpi=150)
