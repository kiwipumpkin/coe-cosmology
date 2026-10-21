import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
S=np.load('mc_sky.npy',mmap_mode='r');sm=np.load('/workspace/redweb2/sm_red.npy');RS=[0.5,0.9,2.0,4.0,8.0,14.0]
t=np.loadtxt('/workspace/joint/lcdm.txt');l=np.arange(2,31);Lt=np.interp(l,t[:,0],t[:,1])
k=851;E=[1,.5,.2,.1,.05,.02,.01,.005,.002,.0005];cols=[plt.cm.turbo(i/8) for i in range(9)]
for R in [4.0,8.0,14.0]:
  iR=RS.index(R);best=min(E,key=lambda e:chi(Dl(sm+e*S[k][iR])))
  plt.figure(figsize=(11,7));plt.plot(l,Lt,'k-',lw=2,label='ΛCDM χ²=20.1');c=0
  for e in [best*f for f in (0.4,0.55,0.7,0.85,1,1.2,1.45,1.75,2.1)]:
    D=Dl(sm+e*S[k][iR]);y=D[2:31];m=np.median(y[2:])
    plt.plot(l,y,'.-',color=cols[c],lw=2.6 if c==4 else 1.4,label=f'{e*100:.3g}%: χ²={chi(D):.0f}, r²={r2(D):.2f}');c+=1
  plt.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='r',ms=5,capsize=2,label='Planck TT')
  plt.ylim(0,2600);plt.xlim(1.5,30.5);plt.xlabel('ℓ');plt.ylabel('D_ℓ (μK²)');plt.title(f'red + Monte Carlo web #851, kernel at {R} Gpc, linear')
  plt.legend(loc='lower center',ncol=3,fontsize=9,framealpha=.85);plt.tight_layout();plt.savefig(f'mc851_R{R}_linfine.png',dpi=100);plt.close()
