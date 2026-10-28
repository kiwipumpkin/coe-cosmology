# Fourier power of connect-the-dots (linear interp, periodic) Migkas A/A_all, dipole removed, m=2..2500, vs Planck 2018 TT l=2..2508.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];o=np.argsort(l);l,y=l[o],y[o]
N=2**16;g=np.arange(N)*360/N;f=np.interp(g,l,y,period=360)
P=np.abs(np.fft.rfft(f)/N)**2*2;M=2500;ms=np.arange(2,M+1);Ps=P[2:M+1];sh=Ps/Ps.sum()
p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');L=p[:,0];Dl=p[:,1];Pl=(2*L+1)*2*np.pi*Dl/(L*(L+1));shp=Pl/Pl.sum()
fig,ax=plt.subplots(2,1,figsize=(11,8))
for a in ax:
  a.plot(ms,sh,'b-',lw=1,label='Migkas connect-the-dots, azimuthal m (dipole removed), share of power')
  a.plot(L,shp,'r-',lw=1,label='Planck 2018 TT, spherical $\\ell$, $(2\\ell+1)C_\\ell$ share')
  a.set_yscale('log');a.set_xlabel('m  /  $\\ell$');a.set_ylabel('fraction of power (m or $\\ell\\geq2$)');a.grid(alpha=.3,which='both');a.set_xlim(2,2508)
ax[0].legend(fontsize=8);ax[0].set_title('Linear m/$\\ell$ axis');ax[1].set_xscale('log');ax[1].set_title('Log m/$\\ell$ axis')
plt.tight_layout();plt.savefig('A_lines_spectrum.png',dpi=140)
print('m2-8',sh[:7].sum(),'9-35',sh[7:34].sum(),'>35',sh[34:].sum(),'m2',sh[0],'m3',sh[1])
