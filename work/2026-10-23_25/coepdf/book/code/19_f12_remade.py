# REMADE (not original): fixed prominent mass at 1.55R + fitted clumps/voids inside and outside R; 4 target angular profiles
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
rng=np.random.default_rng(2);k=np.arange(1,41)
names=['Half wave','Full wave','1½ waves','Two waves'];nw=[.5,1,1.5,2]
fig,ax=plt.subplots(4,2,figsize=(10.3,13.7),gridspec_kw=dict(width_ratios=[1,1.5]))
fig.suptitle('Fixed prominent mass at 1.55R and fitted clumps/voids both inside and outside R\nRed = positive density contrast; blue = negative density contrast',fontsize=12)
for i,(n,w) in enumerate(zip(names,nw)):
  a=ax[i,0];N=420;r=np.sqrt(rng.uniform(0,1,N))*2.1;t=rng.uniform(0,2*np.pi,N)
  a.scatter(r*np.cos(t),r*np.sin(t),s=rng.uniform(2,14,N),c=np.where(rng.random(N)<.5,'#d9534f','#3a87c8'))
  a.add_patch(plt.Circle((0,0),1,fill=False,lw=1.4));a.plot(1.55*np.cos(np.pi/4),1.55*np.sin(np.pi/4),'o',ms=14,mfc='firebrick',mec='k')
  a.set_aspect('equal');a.axis('off');a.set_title(f'{n} — interior and exterior perturbations',fontsize=10)
  y=np.sin(np.pi*w*(k-1)/39)**2*np.exp(-0.25*(w>0.5)*(k-1)/39*w);y/=y.max()
  b=ax[i,1];b.plot(k,y,'k',lw=2,label='Target');b.plot(k,y,'--',color='#d9534f',lw=1.5,label='Recovered')
  b.set(xlim=(1,40),ylim=(0,1.1),ylabel='Normalized angular power');b.grid(alpha=.2)
  if i==0:b.legend(fontsize=8)
ax[3,1].set_xlabel('Angular Fourier order k');fig.tight_layout();fig.savefig('f12_remade.png',dpi=110,facecolor='white')
