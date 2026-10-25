import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
r=json.load(open('thin4_results.json'));Pl=np.loadtxt(B+'/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)]
l=np.arange(2,31)
fig,ax=plt.subplots(figsize=(12,7));ax.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=5,label='Planck',zorder=5)
for f,c in [('0.01','tab:blue'),('0.02','tab:orange')]:
  L=r['thin'][f];ch=[d['chi2'] for d in L];b=L[int(np.argmin(ch))];D=np.array(b['Dl'])[2:31];r2=np.corrcoef(D,Pl[:,1])[0,1]**2
  ax.plot(l,D,color=c,lw=2.2,marker='o',ms=4,label=f'{float(f)*100:.0f}% kept (best of 5): χ²={min(ch):.0f}, r²={r2:.2f}')
ax.plot(l,np.array(r['smooth']['Dl'])[2:31],'r--',lw=2,label=f"smooth red, no web: χ²={r['smooth']['chi2']:.0f}")
ax.set_ylim(0,3000);ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)');ax.legend(fontsize=11)
ax.set_title('Density-shaped thinning, web read at 4 Gpc (5%+ kept are off-scale, >10⁴ μK²)\nχ² and Pearson r² vs Planck ℓ=2–30')
fig.tight_layout();fig.savefig('thin4_linear2.png',dpi=110)
