# Standard D_l = l(l+1)C_l/2pi presentation. Planck 2018 TT D_l (uK^2, left axis, with 68% errors) and
# Migkas A/A_all residual ring harmonics m=2..35 as m(m+1)amp^2/2pi (dimensionless, right axis). Dipole removed.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall']
c=np.loadtxt('/workspace/cosmology-course/notes/catwise-eta-migkas-fit-d14.3.csv',delimiter=',',skiprows=1)
r=y-np.interp(l,c[:,0],c[:,1],period=360);L=np.radians(l);M=(len(l)-1)//2
X=np.stack([np.ones_like(L)]+sum([[np.cos(k*L),np.sin(k*L)] for k in range(1,M+1)],[]),1)
co=np.linalg.lstsq(X,r,rcond=None)[0];a2=co[1::2]**2+co[2::2]**2;m=np.arange(1,M+1);k=m>=2
Dm=m[k]*(m[k]+1)*a2[k]/(2*np.pi)
P=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt')
fig,axs=plt.subplots(2,1,figsize=(11,9))
for ax,xs in zip(axs,['log','linear']):
  ax.errorbar(P[:,0],P[:,1],yerr=[P[:,2],P[:,3]],fmt='-',color='darkorange',ecolor='moccasin',lw=1,label='Planck 2018 TT $D_\\ell$')
  ax.set_xscale(xs);ax.set_ylabel('$D_\\ell=\\ell(\\ell+1)C_\\ell/2\\pi$  [$\\mu$K$^2$]',color='darkorange');ax.set_xlabel('m or $\\ell$ (%s)'%xs);ax.grid(alpha=.3)
  t=ax.twinx();t.plot(m[k],Dm,'o-',color='steelblue',ms=4,label='Migkas residual m(m+1)amp$^2$/2$\\pi$, m=2..35')
  t.set_ylabel('residual m(m+1)amp$^2$/2$\\pi$ (dimensionless)',color='steelblue');t.set_ylim(0,Dm.max()*1.15)
  h1,l1=ax.get_legend_handles_labels();h2,l2=t.get_legend_handles_labels();ax.legend(h1+h2,l1+l2,loc='upper right')
axs[0].set_title('Standard $D_\\ell$ presentation, dipole removed (separate y-axes)')
fig.text(.5,.005,'m = azimuthal harmonic in Galactic longitude along a ring; $\\ell$ = spherical multipole. Not the same quantity. Residual stops at m=35 (Nyquist, 71 cones).',ha='center',fontsize=8)
plt.tight_layout(rect=[0,.02,1,1]);plt.savefig('A_residual_vs_planck_Dl.png',dpi=150)
i=np.argmax(P[:,1]);print(P[i,:2])
