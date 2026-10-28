# Weighted Fourier-series fit to Migkas A/A_all, M=35; dipole removed: shares over m=2..35 vs Planck 2018 TT C_l shares over l=2..35.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];s=.5*(d['A_err_low']+d['A_err_high'])
L=np.radians(l);M=35;X=np.stack([np.ones_like(L)]+sum([[np.cos(k*L),np.sin(k*L)] for k in range(1,M+1)],[]),1)
co=np.linalg.lstsq(X/s[:,None],y/s,rcond=None)[0];a2=(co[1::2]**2+co[2::2]**2)[1:];m=np.arange(2,M+1);f=a2/a2.sum()
P=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt')[:34];C=P[:,1]*2*np.pi/(P[:,0]*(P[:,0]+1));g=C/C.sum()
fig,ax=plt.subplots(figsize=(9,5.5));ax.loglog(m,f,'o-',color='royalblue',label='Migkas Fourier fit M=35, amp$^2$ share (m=2..35)')
ax.loglog(P[:,0],g,'s-',color='darkorange',label='Planck 2018 TT, $C_\\ell$ share ($\\ell$=2..35)')
ax.set(xlabel='m or $\\ell$',ylabel='fraction of power over 2..35 (log)',title='Dipole removed');ax.grid(alpha=.3,which='both');ax.legend(fontsize=8)
fig.text(.5,.005,'m = azimuthal ring harmonic; $\\ell$ = spherical multipole. Not the same quantity. Weighted by 1/sigma^2.',ha='center',fontsize=8)
plt.tight_layout(rect=[0,.03,1,1]);plt.savefig('A_fourier_fit_power_noDipole_M35.png',dpi=140)
print(np.round(100*f[:7],1),np.round(100*g[:7],1),round(100*f[7:].sum(),1),round(100*g[7:].sum(),1))
