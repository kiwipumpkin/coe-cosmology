# Power spectrum of weighted Fourier-series fit to Migkas A/A_all (M=35), dipole KEPT; fractions over m=1..8.
# Planck 2018 TT l=2..8 from file plus l=1 dipole A=3362.08 uK (Planck 2018 I, arXiv:1807.06205): C1=4*pi*A^2/9. Fractions over l=1..35.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];s=.5*(d['A_err_low']+d['A_err_high'])
L=np.radians(l);M=35;X=np.stack([np.ones_like(L)]+sum([[np.cos(k*L),np.sin(k*L)] for k in range(1,M+1)],[]),1)
co=np.linalg.lstsq(X/s[:,None],y/s,rcond=None)[0];a2=co[1::2]**2+co[2::2]**2;m=np.arange(1,M+1);f=a2/a2.sum()
P=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt')[:34];ll=np.r_[1,P[:,0]];C=np.r_[4*np.pi*3362.08**2/9,P[:,1]*2*np.pi/(P[:,0]*(P[:,0]+1))];g=C/C.sum()
fig,ax=plt.subplots(figsize=(9,5.5));ax.semilogy(m,f,'o-',color='royalblue',label='Migkas Fourier fit M=35, amp$^2$ share (m=1..35)')
ax.semilogy(ll,g,'s-',color='darkorange',label='Planck 2018 TT incl. dipole, $C_\\ell$ share ($\\ell$=1..35)')
ax.set(xlabel='m or $\\ell$',ylabel='fraction of power over 1..35 (log)',title='Dipole kept');ax.grid(alpha=.3,which='both');ax.legend(fontsize=8);ax.set_xscale('log')
fig.text(.5,.005,'m = azimuthal ring harmonic; $\\ell$ = spherical multipole. Not the same quantity. Planck dipole is largely kinematic.',ha='center',fontsize=8)
plt.tight_layout(rect=[0,.03,1,1]);plt.savefig('A_fourier_fit_power_withDipole_M35.png',dpi=140)
print(np.round(100*f,2));print(np.array2string(100*g,precision=5))
