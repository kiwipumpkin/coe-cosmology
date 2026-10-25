# REMADE (not original): non-axisymmetric source shell r=1.012R, Gaussian peak at l=220 (sigma 32, 8)
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
rng=np.random.default_rng(1);ell=np.arange(100,331)
fig=plt.figure(figsize=(12,10));fig.suptitle('Non-axisymmetric source shell: density varies across latitude AND longitude\nExterior shell r=1.012R • gravitational potential, not CMB temperature',fontsize=13)
for k,s in enumerate((32,8)):
  ax=fig.add_subplot(2,2,2*k+1,projection='3d');u,v=np.mgrid[0:np.pi:120j,0:2*np.pi:240j]
  f=rng.standard_normal(u.shape);clip=np.percentile(abs(f),97);col=plt.cm.RdBu_r((np.clip(f,-clip,clip)/clip+1)/2)
  ax.plot_surface(np.sin(u)*np.cos(v),np.sin(u)*np.sin(v),np.cos(u),facecolors=col,rstride=1,cstride=1,linewidth=0,shade=False)
  ax.set_box_aspect((1,1,1));ax.axis('off');ax.set_title(f'Distributed 3D clumps and voids\nσℓ={s}; 97th-percentile color clip')
  a=fig.add_subplot(2,2,2*k+2);t=np.exp(-.5*((ell-220)/s)**2)
  a.plot(ell,t,'k',lw=2.2,label='Prescribed target');a.plot(ell,t,'--',color='#d24a3b',lw=1.8,label='Calculated spherical power')
  a.axvline(220,color='.4',ls=':',lw=1);a.set(xlim=(110,330),ylim=(0,1.08),xlabel='Spherical multipole ℓ',ylabel='Normalized potential Dℓ',title=f'Peak ℓ=220, σℓ={s}');a.grid(alpha=.2);a.legend(fontsize=8)
fig.tight_layout();fig.savefig('f13_remade.png',dpi=110,facecolor='white')
