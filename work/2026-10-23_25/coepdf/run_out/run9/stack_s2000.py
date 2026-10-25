import json,os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
os.chdir('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained')
Pl=np.loadtxt('../../../reproduce_cmb/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)]
th=np.loadtxt('/tmp/th.txt');th=th[(th[:,0]>=2)&(th[:,0]<=30)]
S=(1,.5,.25,.1,.05,.02,.01,.005,.002,.001)
d=json.load(open('red_web_s2000_results.json'))
rows=[('0.23','red_surface'),('0.6','R0.6'),('0.8','R0.8'),('0.9','R0.9'),('1','R1'),('2','R2'),('3','R3'),('4','R4'),('14.3','paper_R14.3')]
fig,ax=plt.subplots(len(rows),1,figsize=(10,3.6*len(rows)),sharex=True);l=np.arange(31);cm=plt.cm.viridis
L=Pl[:,0].astype(int);y=Pl[:,1];r2=lambda p:1-np.sum((y-p)**2)/np.sum((y-y.mean())**2);best=[]
for a,(name,tag) in zip(ax,rows):
  a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=3,label='Planck')
  a.plot(th[:,0],th[:,1],'k-',lw=1.5,label='ΛCDM χ²=20.1')
  sm=d['smooth'];a.plot(l[2:],np.array(sm['Dl'])[2:],'r--',lw=1.5,label=f"smooth red χ²={sm['chi2']:.1f}")
  for i,s in enumerate(S):
    r=d[f'{tag}_{s}'];D=np.array(r['Dl']);best.append((r['chi2'],r2(D[L]),name,s));a.plot(l[2:],np.maximum(D[2:],1e-3),color=cm(i/9),lw=1,label=f'{s*100:g}% χ²={r["chi2"]:.4g}')
  a.set_yscale('log');a.set_ylim(1e-3,1e7);a.set_title(f'red + web (seed 2000), kernel at {name} Gpc');a.set_ylabel('D_ℓ (μK²)');a.legend(fontsize=6.5,ncol=3,loc='upper right')
ax[-1].set_xlim(1.5,30.5);ax[-1].set_xlabel('ℓ');fig.tight_layout();fig.savefig('red_web_s2000_stack.png',dpi=110)
best.sort();[print(b) for b in best[:6]]
D=np.array(d['R0.8_0.05']['Dl']);plt.figure(figsize=(10,5));plt.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=4,label='Planck')
plt.plot(l[2:],D[2:],'-o',color='teal',ms=3,label=f"red + web seed 2000, 0.8 Gpc, 5% (χ²={d['R0.8_0.05']['chi2']:.1f}, R²={r2(D[L]):.2f})")
plt.xlabel('ℓ');plt.ylabel('D_ℓ (μK²)');plt.xlim(1.5,30.5);plt.ylim(0);plt.legend();plt.title('Linear: Planck vs red + web seed 2000 (0.8 Gpc, 5%)');plt.tight_layout();plt.savefig('red_web_s2000_08_5pct_linear.png',dpi=120)
