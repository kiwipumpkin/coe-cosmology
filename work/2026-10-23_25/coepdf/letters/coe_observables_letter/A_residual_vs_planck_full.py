# Migkas A/A_all residual ring harmonics m=2..35 (Nyquist for 71 cones at 5 deg) vs Planck 2018 TT C_l, l=2..lmax. Dipole (m=1, l=1) removed.
# Each spectrum normalized to sum 1 over its own range (2..35 for m, 2..lmax for l).
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall']
c=np.loadtxt('/workspace/cosmology-course/notes/catwise-eta-migkas-fit-d14.3.csv',delimiter=',',skiprows=1)
r=y-np.interp(l,c[:,0],c[:,1],period=360);L=np.radians(l);N=len(l);M=(N-1)//2
X=np.stack([np.ones_like(L)]+sum([[np.cos(k*L),np.sin(k*L)] for k in range(1,M+1)],[]),1)
coef=np.linalg.lstsq(X,r,rcond=None)[0];amp2=coef[1::2]**2+coef[2::2]**2
m=np.arange(1,M+1);k=m>=2;pm=amp2[k]/amp2[k].sum()
P=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');ll=P[:,0];C=P[:,1]*2*np.pi/(ll*(ll+1));pl=C/C.sum()
fig,ax=plt.subplots(figsize=(11,5.5))
ax.loglog(ll,pl,'-',color='darkorange',lw=1,label='Planck 2018 TT $C_\\ell$, $\\ell$=2..%d'%ll.max())
ax.loglog(m[k],pm,'o-',color='steelblue',label='Migkas $A/A_{\\rm all}$ residual, ring harmonic m=2..%d (Nyquist, %d cones)'%(M,N))
ax.set(xlabel='m or $\\ell$',ylabel='fractional power (each sums to 1 over its own range)',title='Dipole removed: residual ring harmonics vs Planck TT')
ax.grid(alpha=.3,which='both');ax.legend()
fig.text(.5,.005,'m = azimuthal harmonic in Galactic longitude along a ring; $\\ell$ = spherical multipole. Not the same quantity.',ha='center',fontsize=8)
plt.tight_layout(rect=[0,.03,1,1]);plt.savefig('A_residual_vs_planck_full.png',dpi=150)
print('N',N,'M',M,'lmax',ll.max(),'rank',np.linalg.matrix_rank(X));print(np.round(100*pm,2))
