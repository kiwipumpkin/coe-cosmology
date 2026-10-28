# Each truncated connect-the-dots Fourier sum (m<=2,8,35,2000): fractional power per m (m>=1, dipole kept) overlaid with Planck 2018 TT (2l+1)C_l share over the same range (l=1 from dipole A=3362.08 uK).
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];o=np.argsort(l);l,y=l[o],y[o]
N=2**16;g=np.arange(N)*360/N;f=np.interp(g,l,y,period=360);A2=4*np.abs(np.fft.rfft(f)/N)**2
p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');L=np.r_[1,p[:,0]];Ad=3362.08;C1=4*np.pi/9*Ad**2
Cl=np.r_[C1,2*np.pi*p[:,1]/(p[:,0]*(p[:,0]+1))];Pl=(2*L+1)*Cl
fig,axs=plt.subplots(2,2,figsize=(12,8))
for a,(M,cl) in zip(axs.flat,[(2,'orange'),(8,'green'),(35,'purple'),(2000,'b')]):
  m=np.arange(1,M+1);pm=A2[1:M+1];pm=pm/pm.sum();k=L<=M;pp=Pl[k]/Pl[k].sum()
  a.loglog(m,pm,'o-' if M<=35 else '-',color=cl,ms=4,lw=1,label=f'Migkas lines, m=1..{M}')
  a.loglog(L[k],pp,'s-' if M<=35 else '-',color='red',ms=4,lw=1,label=f'Planck TT, $\\ell$=1..{M}')
  a.set_title(f'm, $\\ell$ = 1..{M} (each normalised to 1)');a.set_xlabel('m or $\\ell$');a.set_ylabel('fraction of power');a.grid(alpha=.3,which='both');a.set_xlim(.8,max(M*1.3,3));a.legend(fontsize=8)
fig.text(.5,.005,'m = azimuthal ring harmonic; $\\ell$ = spherical multipole. Not the same quantity. Planck dipole is largely kinematic.',ha='center',fontsize=9)
plt.tight_layout(rect=(0,.03,1,1));plt.savefig('A_lines_harm_spectra_planck.png',dpi=140)
