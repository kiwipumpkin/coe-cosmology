# Fourier power of the periodic cubic spline through Migkas A/A_all (dipole m=1 removed), m=2..2500, vs Planck 2018 TT l=2..2508.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall']
o=np.argsort(l);l,y=l[o],y[o]
s=CubicSpline(np.r_[l,l[0]+360],np.r_[y,y[0]],bc_type='periodic')
N=2**16;g=np.arange(N)*360/N;f=s(np.where(g<l[0],g+360,g))
c=np.fft.rfft(f)/N;P=np.abs(c)**2*2;m=np.arange(len(P))
M=2500;ms=m[2:M+1];Ps=P[2:M+1];sh=Ps/Ps.sum()
p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');L=p[:,0];Dl=p[:,1];Cl=2*np.pi*Dl/(L*(L+1));Pl=(2*L+1)*Cl;shp=Pl/Pl.sum()
fig,ax=plt.subplots(2,1,figsize=(11,8))
for a in ax:
  a.plot(ms,sh,'b-',lw=1,label='Migkas spline, azimuthal m (dipole removed), share of power')
  a.plot(L,shp,'r-',lw=1,label='Planck 2018 TT, spherical $\\ell$, $(2\\ell+1)C_\\ell$ share')
  a.set_yscale('log');a.set_xlabel('m  /  $\\ell$');a.set_ylabel('fraction of power (m or $\\ell\\geq2$)');a.grid(alpha=.3,which='both')
ax[0].set_xlim(2,2508);ax[0].legend(fontsize=8);ax[0].set_title('Linear m/$\\ell$ axis')
ax[1].set_xscale('log');ax[1].set_xlim(2,2508);ax[1].set_title('Log m/$\\ell$ axis')
ax[0].set_ylim(1e-14,1)
plt.tight_layout();plt.savefig('A_spline_spectrum.png',dpi=140)
for k in [2,3,4,8,35,72,144,500,1000,2500]:print(k,sh[k-2])
print('share m2-8',sh[:7].sum(),'9-35',sh[7:34].sum(),'>35',sh[34:].sum())
