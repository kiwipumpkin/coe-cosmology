import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
lon=np.linspace(-np.pi,np.pi,360); lat=np.linspace(-np.pi/2,np.pi/2,180)
L,B=np.meshgrid(lon,lat)
# Y20-like + Y22-like real quadrupole
Q=0.5*(3*np.sin(B)**2-1)+0.6*np.cos(B)**2*np.cos(2*L)
fig=plt.figure(figsize=(16,9))
for i,(f,t) in enumerate([(Q,'Q'),(-Q,'−Q (exact opposite)'),(Q-Q,'Q + (−Q) = 0')]):
    ax=fig.add_subplot(2,3,i+1,projection='mollweide')
    im=ax.pcolormesh(L,B,f,cmap='RdBu_r',vmin=-1,vmax=1,shading='auto')
    ax.set_title(t,fontsize=18); ax.set_xticklabels([]); ax.set_yticklabels([])
ax=fig.add_subplot(2,1,2)
th=np.degrees(lon); q=Q[:,0]; 
ql=0.5*(3*np.sin(lon/2)**2-1)  # placeholder unused
cut=0.5*(3*np.sin(np.radians(np.linspace(-90,90,181)))**2-1)+0.6*np.cos(np.radians(np.linspace(-90,90,181)))**2
x=np.linspace(-90,90,181)
ax.plot(x,cut,'r',lw=3,label='Q along meridian l=0°')
ax.plot(x,-cut,'b',lw=3,label='−Q along meridian')
ax.plot(x,cut-cut,'k--',lw=3,label='sum = 0')
ax.set_xlabel('Galactic latitude b [deg]',fontsize=16); ax.set_ylabel('ΔT (arb.)',fontsize=16); ax.tick_params(labelsize=14)
ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.18),ncol=3,fontsize=14)
fig.text(0.5,0.5,r'$Q(\hat n)=\sum_m a_{2m}Y_{2m}(\hat n)$,   $-Q$ has $a_{2m}\to-a_{2m}$,   $C_2=\frac{1}{5}\sum_m|a_{2m}|^2$ of sum $=0$',ha='center',fontsize=16)
plt.tight_layout(rect=[0,0,1,0.97]); plt.savefig('opposite_quadrupoles.png',dpi=110)
