# l=1 / m=1 removed (CMB dipole treated as our peculiar motion). Fractional power renormalized over 2..8.
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
M=json.load(open('A_fit_residual_fourier.json'))['modes'];m=np.array([x['m'] for x in M]);a2=np.array([x['amp']**2 for x in M])
k=m>=2;m=m[k];r=a2[k]/a2[k].sum()
d=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');d=d[(d[:,0]>=2)&(d[:,0]<=8)]
l=d[:,0];f=2*np.pi/(l*(l+1));C=d[:,1]*f;S=C.sum();p=C/S;lo=d[:,2]*f/S;hi=d[:,3]*f/S
fig,ax=plt.subplots(1,2,figsize=(12,4.5));w=.38
ax[0].bar(m-w/2,100*r,w,label='A/A_all residual, m (ring harmonic)',color='steelblue')
ax[0].bar(l+w/2,100*p,w,yerr=[100*lo,100*hi],capsize=3,label='Planck 2018 TT, ℓ (C_ℓ)',color='darkorange')
ax[0].set(xlabel='m or ℓ',ylabel='% of power over 2..8',title='Dipole removed: shares over 2..8');ax[0].legend()
ax[1].semilogy(m,r,'o-',label='residual m');ax[1].errorbar(l,p,yerr=[np.minimum(lo,p*.99),hi],fmt='s-',capsize=3,label='Planck ℓ')
ax[1].set(xlabel='m or ℓ',ylabel='fractional power (log)',title='Log comparison');ax[1].legend()
for x in ax:x.set_xticks(range(2,9));x.grid(alpha=.3)
fig.text(.5,.005,'m = azimuthal harmonic in Galactic longitude along a ring; ℓ = spherical multipole. Not the same quantity.',ha='center',fontsize=8)
plt.tight_layout(rect=[0,.03,1,1]);plt.savefig('A_residual_power_vs_planck_noDipole.png',dpi=150)
for i in range(7):print(int(m[i]),round(100*r[i],1),int(l[i]),round(100*p[i],1))
