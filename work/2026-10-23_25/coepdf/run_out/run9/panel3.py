import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
Pl=np.loadtxt(B+'/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)];l=np.arange(2,31)
fo=json.load(open('force4_results.json'));t=json.load(open('thin4_results.json'))
def r2(D):return np.corrcoef(D,Pl[:,1])[0,1]**2
sets=[]
sets.append(('Original web, all bodies, strength sweep (4 Gpc)',[(f'strength {float(k)*100:g}%',v['Dl'],v['chi2']) for k,v in sorted(fo['orig'].items(),key=lambda x:float(x[0]))]))
sets.append(('Density-following web (seed 1), all bodies, strength sweep (4 Gpc)',[(f'strength {float(k)*100:g}%',v['Dl'],v['chi2']) for k,v in sorted(fo['densweb_s1'].items(),key=lambda x:float(x[0]))]))
th=[]
for k,L in sorted(t['thin'].items(),key=lambda x:float(x[0])):
  b=min(L,key=lambda d:d['chi2']);th.append((f'{float(k)*100:g}% kept (best of {len(L)})',b['Dl'],b['chi2']))
sets.append(('Original web, density-shaped thinning, remove bodies (4 Gpc)',th))
fig,axs=plt.subplots(3,1,figsize=(11,27))
for ax,(title,cs) in zip(axs,sets):
  cm=plt.cm.turbo
  for i,(lab,D,chi) in enumerate(cs):
    D=np.array(D)[2:31];ax.plot(l,D,color=cm(i/max(1,len(cs)-1)),lw=1.6,marker='.',label=f'{lab}: χ²={chi:.3g}, r²={r2(D):.2f}')
  ax.plot(l,np.array(t['smooth']['Dl'])[2:31],'k--',lw=1.6,label=f"smooth red, no web: χ²={t['smooth']['chi2']:.1f}")
  ax.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='red',ms=5,label='Planck',zorder=10)
  ax.set_yscale('log');ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)');ax.set_title(title)
  ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.08),ncol=2,fontsize=9,frameon=False)
fig.suptitle('χ² and Pearson r² vs Planck ℓ=2–30',y=0.995)
fig.tight_layout(h_pad=4);fig.savefig('panel3_4gpc.png',dpi=90,bbox_inches='tight')
