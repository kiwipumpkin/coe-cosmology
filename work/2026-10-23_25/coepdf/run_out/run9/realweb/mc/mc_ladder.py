import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
S=np.load('mc_sky.npy',mmap_mode='r');sm=np.load('/workspace/redweb2/sm_red.npy');RS=[0.5,0.9,2.0,4.0,8.0,14.0]
t=np.loadtxt('/workspace/joint/lcdm.txt');l=np.arange(2,31);Lt=np.interp(l,t[:,0],t[:,1])
k=int(np.random.default_rng(2026).integers(len(S)));E=[1,.5,.2,.1,.05,.02,.01,.005,.002,.0005]
fig,ax=plt.subplots(len(RS),1,figsize=(11,5*len(RS)))
for a,iR,R in zip(ax,range(6),RS):
  a.plot(l,Lt,'k-',lw=2,label='ΛCDM χ²=20.1')
  D0=Dl(sm);a.plot(l,D0[2:31],'r--',lw=2,label=f'smooth red χ²={chi(D0):.1f}')
  for j,e in enumerate(E):
    D=Dl(sm+e*S[k][iR]);a.plot(l,D[2:31],color=plt.cm.viridis(j/(len(E)-1)),label=f'{e*100:g}% χ²={chi(D):.4g}')
  a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=4,label='Planck')
  a.set_yscale('log');a.set_ylim(1e-3,None);a.set_ylabel('D_ℓ (μK²)');a.set_title(f'red + Monte Carlo web #{k}, kernel at {R} Gpc')
  a.legend(loc='upper right',ncol=3,fontsize=7)
ax[-1].set_xlabel('ℓ');plt.tight_layout();plt.savefig('mc_ladder.png',dpi=90);print(k)
