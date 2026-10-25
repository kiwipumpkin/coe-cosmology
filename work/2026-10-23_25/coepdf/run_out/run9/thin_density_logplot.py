import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
r=json.load(open('thin_density_results.json'));Pl=np.loadtxt(B+'/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)]
l=np.arange(2,31);rows=[]
fig,ax=plt.subplots(figsize=(13,8));ax.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=4,label='Planck',zorder=5)
cm=plt.cm.viridis
for i,(f,L) in enumerate(r['thin'].items()):
  c=[d['chi2'] for d in L];b=L[int(np.argmin(c))];D=np.array(b['Dl'])[2:31];r2=np.corrcoef(D,Pl[:,1])[0,1]**2
  rows.append((f,np.mean(c),min(c),max(c),np.std(c),r2,D[0],len(L)))
  ax.plot(l,D,color=cm(i/12),lw=1.4,marker='.',label=f'f={f}: χ²={min(c):.1f}, r²={r2:.2f}')
ax.plot(l,np.array(r['smooth']['Dl'])[2:31],'r--',label=f"smooth red only χ²={r['smooth']['chi2']:.1f}")
ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)');ax.set_yscale('log');ax.legend(fontsize=8,ncol=2)
ax.set_title('Red universe + seed-2000 web, density-modulated by red ρ, thinned to kept fraction f (full mass), web read at 0.9 Gpc; best of draws per f\nr² = Pearson r² vs Planck ℓ=2–30')
fig.tight_layout();fig.savefig('thin_density_log.png',dpi=110)
for x in rows:print('%s mean %.1f best %.1f worst %.1f sd %.1f r2 %.3f D2 %.0f n=%d'%x)
print({k:v['chi2'] for k,v in r.items() if k.startswith('scaled')})
