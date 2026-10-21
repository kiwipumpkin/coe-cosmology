import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
R=json.load(open('thin4_jump.json'))['rows'];Pl=np.loadtxt(B+'/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)];l=np.arange(2,31)
fig,ax=plt.subplots(figsize=(12,7));ax.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=5,label='Planck',zorder=10)
cols=['tab:gray','tab:red','tab:blue','tab:green','tab:purple','tab:orange'];ls=['-','-','--','-.',':','-']
for i,f in enumerate([0.017,0.0171,0.0172,0.0173,0.0174,0.0175]):
  z=min(R,key=lambda q:abs(q['f']-f));D=np.array(z['Dl'])[2:31]
  ax.plot(l,D,color=cols[i],ls=ls[i],lw=2.4-0.25*i,marker='o',ms=3,label=f"{f*100:.2f}% (+{z['n_added']} bodies): χ²={z['chi2']:.1f}, r²={z['r2']:.3f}")
  print(f,z['n_added'],round(z['chi2'],2),round(z['r2'],4),round(D[-1]))
ax.set_ylim(0,3000);ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)');ax.legend(fontsize=10)
ax.set_title('Draw 0, density-shaped thinning, 4 Gpc read: steps 1.70%→1.75%\nχ² and Pearson r² vs Planck ℓ=2–30')
fig.tight_layout();fig.savefig('thin4_jump_steps.png',dpi=110)
