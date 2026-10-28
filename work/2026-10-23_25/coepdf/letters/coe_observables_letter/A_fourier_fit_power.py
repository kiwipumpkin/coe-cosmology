# Power spectrum of the weighted Fourier-series fit to Migkas A/A_all vs l (sigma = mean asym error), for M=2,3,8.
# Power per mode = amp^2 (and D-style m(m+1)amp^2/2pi), dipole m=1 shown but excluded from the fractional normalisation; vs Planck 2018 TT.
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];s=.5*(d['A_err_low']+d['A_err_high'])
L=np.radians(l);P=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');ll=P[:,0];C=P[:,1]*2*np.pi/(ll*(ll+1))
def fit(M):
  X=np.stack([np.ones_like(L)]+sum([[np.cos(k*L),np.sin(k*L)] for k in range(1,M+1)],[]),1)
  co=np.linalg.lstsq(X/s[:,None],y/s,rcond=None)[0];return co[1::2]**2+co[2::2]**2
out={};fig,ax=plt.subplots(1,3,figsize=(16,5))
for M,c in [(2,'crimson'),(3,'k'),(8,'royalblue')]:
  a2=fit(M);m=np.arange(1,M+1);out[M]={'amp2':a2.tolist()}
  ax[0].semilogy(m,a2,'o-',color=c,label='M=%d'%M)
  ax[1].plot(m,m*(m+1)*a2/(2*np.pi),'o-',color=c,label='M=%d'%M)
  if M>=2:k=m>=2;f=a2[k]/a2[k].sum();out[M]['frac_m2up']=f.tolist();ax[2].loglog(m[k],f,'o-',color=c,label='Migkas fit M=%d'%M)
ax[2].loglog(ll[:7],C[:7]/C[:7].sum(),'s-',color='darkorange',label='Planck $C_\\ell$ (frac over 2..8)')
ax[0].set(xlabel='m',ylabel='amp$^2$',title='Fit power per mode (incl. dipole m=1)');ax[1].set(xlabel='m',ylabel='m(m+1)amp$^2$/2$\\pi$',title='D-style')
ax[2].set(xlabel='m or $\\ell$',ylabel='fractional power, m$\\geq$2',title='Dipole removed vs Planck')
for a in ax:a.grid(alpha=.3,which='both');a.legend(fontsize=8)
fig.text(.5,.005,'m = azimuthal harmonic in Galactic longitude along a ring; $\\ell$ = spherical multipole. Not the same quantity. Weighted by 1/sigma^2.',ha='center',fontsize=8)
plt.tight_layout(rect=[0,.03,1,1]);plt.savefig('A_fourier_fit_power.png',dpi=140);json.dump(out,open('A_fourier_fit_power.json','w'),indent=1)
for M in out:print(M,np.round(out[M]['amp2'],6),np.round(100*np.array(out[M]['frac_m2up']),1))
