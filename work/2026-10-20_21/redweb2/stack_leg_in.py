import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
Pl=np.loadtxt(B+'/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)];l=np.arange(2,31)
fo=json.load(open('force4_results.json'));t=json.load(open('thin4_results.json'))
r2=lambda D:np.corrcoef(D,Pl[:,1])[0,1]**2
srt=lambda d:sorted(d.items(),key=lambda x:float(x[0]))
P=[('seed-2000 web, all bodies, strength sweep, 4 Gpc',[(f'strength {k}',v['Dl'],v['chi2']) for k,v in srt(fo['orig'])]),
   ('density-following web seed 1, strength sweep, 4 Gpc',[(f'strength {k}',v['Dl'],v['chi2']) for k,v in srt(fo['densweb_s1'])])]
th=[]
for k,L in srt(t['thin']):
  b=min(L,key=lambda d:d['chi2']);th.append((f'kept {k} (best of {len(L)})',b['Dl'],b['chi2']))
for k in ['0.01','0.0175','0.02']:
  v=t['scaled_'+k];th.append((f'strength {k} (all density-kept bodies)',v['Dl'],v['chi2']))
P.append(('Density-shaped thinning, web read at 4 Gpc',th))
fig,axs=plt.subplots(3,1,figsize=(12,30))
for ax,(ti,cs) in zip(axs,P):
  for i,(lab,D,c) in enumerate(cs):
    D=np.array(D)[2:31];ax.plot(l,D,color=plt.cm.viridis(i/max(1,len(cs)-1)),lw=1.2,label=f'{lab}: χ²={c:.4g}, r²={r2(D):.2f}')
  S=np.array(t['smooth']['Dl'])[2:31];ax.plot(l,S,'k-',lw=2.5,label=f"smooth red (no web): χ²={t['smooth']['chi2']:.1f}, r²={r2(S):.2f}")
  ax.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='red',ms=5,label='Planck TT',zorder=10)
  ax.set_yscale('log');y0,y1=ax.get_ylim();ax.set_ylim(1e-3,y1*2)
  ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)');ax.set_title(ti)
  ax.legend(loc='lower center',ncol=2,fontsize=8,framealpha=0.9)
fig.tight_layout(h_pad=3);fig.savefig('stack_orig3_legin.png',dpi=90)
