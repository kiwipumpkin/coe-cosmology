# SIN-AMP variant: residuals from A_fit_residual_fourier_sinamp.json (model*(1-Amp sin(l+phi)), Amp=-0.06, phi=0.35).
# Residual harmonic power (m=1..8) vs Planck 2018 CMB TT low multipoles.
# Planck TT: COM_PowerSpect_CMB-TT-full_R3.01.txt (IRSA mirror of Planck Legacy Archive)
#  https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/COM_PowerSpect_CMB-TT-full_R3.01.txt
# Dipole: 3362.08 +/- 0.99 uK, Planck 2018 results I, arXiv:1807.06205.
# Pure dipole T=A cos(theta): C_1 = 4*pi*A^2/9, D_1 = l(l+1)C_1/2pi = 4A^2/9.
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=json.load(open('A_fit_residual_fourier_sinamp.json'));M=d['modes']
m=np.array([x['m'] for x in M]);P=np.array([x['amp']**2 for x in M]);fr=P/P.sum()
t=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');t=t[(t[:,0]>=2)&(t[:,0]<=8)]
l=t[:,0];Dl=t[:,1];lo=t[:,2];hi=t[:,3];fac=2*np.pi/(l*(l+1))
Cl=Dl*fac;Clo=lo*fac;Chi=hi*fac;S=Cl.sum();fC=Cl/S
A1=3362.08;D1=4*A1**2/9;C1=np.pi*D1
l18=np.r_[1,l];C18=np.r_[C1,Cl];f18=C18/C18.sum()
fig,ax=plt.subplots(2,2,figsize=(12,8.5));a=ax.ravel()
a[0].bar(m,fr,color='steelblue');[a[0].text(x,y,'%.1f%%'%(100*y),ha='center',va='bottom',fontsize=8) for x,y in zip(m,fr)]
a[0].set(xlabel='azimuthal harmonic m (Galactic longitude, along ring)',ylabel='fraction of power',title='A/A_all residual (sin-amp): amp$^2$ / sum, m=1..8')
a[1].bar(l,fC,color='darkorange',yerr=[Clo/S,Chi/S],capsize=3);[a[1].text(x,y,'%.1f%%'%(100*y),ha='center',va='bottom',fontsize=8) for x,y in zip(l,fC)]
a[1].set(xlabel='spherical multipole l',ylabel='fraction of C_l summed over l=2..8',title='Planck 2018 TT: C_l / sum, l=2..8 (68% errors)')
a[2].semilogy(l18,f18,'s-',color='darkred');a[2].set(xlabel='spherical multipole l',ylabel='C_l / sum over l=1..8 (log)',title='Planck 2018 with dipole (l=1: A=3362 uK)')
a[3].semilogy(m,fr,'o-',color='steelblue',label='residual, m=1..8 (sum=1)');a[3].semilogy(l,fC,'s-',color='darkorange',label='Planck C_l, l=2..8 (sum=1)')
a[3].set(xlabel='m or l',ylabel='normalized power (log)',title='Comparison of normalized spectra');a[3].legend(fontsize=8)
for x in a:x.set_xticks(range(1,9));x.grid(alpha=.3)
fig.text(.5,.005,'Note: m is an azimuthal harmonic in Galactic longitude along a ring; l is a spherical multipole. They are not the same quantity.',ha='center',fontsize=9,style='italic')
plt.tight_layout(rect=[0,.02,1,1]);plt.savefig('A_residual_power_vs_planck_sinamp.png',dpi=150)
for i in range(len(l)):print(int(l[i]),Dl[i],lo[i],hi[i],'%.3e'%Cl[i],'%.4f'%fC[i])
print('dipole D1=%.4e C1=%.4e frac(l=1..8)=%.6f'%(D1,C1,f18[0]))
