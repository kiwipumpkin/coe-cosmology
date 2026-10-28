# Treat Migkas A/A_all vs l as a smooth curve; bin error-bar top (y+eh) and bottom (y-el) taken as the smooth system's max and min.
# Fourier (ring harmonics m, up to Nyquist 35 for 71 cones) of upper envelope, lower envelope, and central values; dipole m=1 removed.
# Shown in D-style m(m+1)amp^2/2pi (right axis) against Planck 2018 TT D_l (left axis), plus fractional power on log-log.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high']
L=np.radians(l);M=(len(l)-1)//2;m=np.arange(1,M+1)
X=np.stack([np.ones_like(L)]+sum([[np.cos(k*L),np.sin(k*L)] for k in range(1,M+1)],[]),1)
def spec(v):co=np.linalg.lstsq(X,v,rcond=None)[0];return co[1::2]**2+co[2::2]**2
S={'upper envelope (max)':(spec(y+eh),'crimson'),'central values':(spec(y),'k'),'lower envelope (min)':(spec(y-el),'royalblue')}
P=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');ll=P[:,0];C=P[:,1]*2*np.pi/(ll*(ll+1))
k=m>=2;fig,ax=plt.subplots(3,1,figsize=(11,13))
for a,xs in zip(ax[:2],['log','linear']):
  a.plot(ll,P[:,1],color='darkorange',lw=1,label='Planck 2018 TT $D_\\ell$');a.set_xscale(xs);a.grid(alpha=.3)
  a.set_ylabel('$D_\\ell$ [$\\mu$K$^2$]',color='darkorange');a.set_xlabel('m or $\\ell$ (%s)'%xs);t=a.twinx();mx=0
  for n,(s,c) in S.items():D=m[k]*(m[k]+1)*s[k]/(2*np.pi);mx=max(mx,D.max());t.plot(m[k],D,'o-',color=c,ms=3,lw=1,label='Migkas '+n)
  t.set_ylim(0,mx*1.1);t.set_ylabel('Migkas m(m+1)amp$^2$/2$\\pi$');h1,l1=a.get_legend_handles_labels();h2,l2=t.get_legend_handles_labels();a.legend(h1+h2,l1+l2,fontsize=8,loc='upper right')
ax[0].set_title('Smooth-envelope assumption: Fourier of Migkas $A/A_{\\rm all}$ max/min envelopes vs Planck (dipole removed)')
ax[2].loglog(ll,C/C.sum(),color='darkorange',lw=1,label='Planck $C_\\ell$')
for n,(s,c) in S.items():ax[2].loglog(m[k],s[k]/s[k].sum(),'o-',color=c,ms=3,lw=1,label='Migkas '+n)
ax[2].set(xlabel='m or $\\ell$',ylabel='fractional power over m,$\\ell\\geq2$');ax[2].grid(alpha=.3,which='both');ax[2].legend(fontsize=8)
fig.text(.5,.004,'m = azimuthal harmonic in Galactic longitude along a ring; $\\ell$ = spherical multipole. Not the same quantity. Migkas stops at m=35 (Nyquist, 71 cones).',ha='center',fontsize=8)
plt.tight_layout(rect=[0,.015,1,1]);plt.savefig('A_envelope_vs_planck.png',dpi=140)
for n,(s,c) in S.items():p=s[k]/s[k].sum();print(n,'m2..8 %:',np.round(100*p[:7],1))
