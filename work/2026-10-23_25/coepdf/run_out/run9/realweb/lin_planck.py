import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb';Pl=np.loadtxt(B+'/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)]
D=json.load(open('realweb_results.json'))['R'];l=np.arange(2,31)
cols=plt.cm.tab10.colors;ls={'a':'--','b':'-','c':':'};nm={'a':'remove','b':'force','c':'mix'}
for R in ['2.0','4.0','8.0','14.0']:
  plt.figure(figsize=(11,7));k=0
  for m in 'abc':
    for e,v in D[R][m].items():
      r=v if isinstance(v,dict) else min(v,key=lambda x:x['chi2'])
      y=np.array(r['Dl'])[2:31];med=np.median(y[2:])
      if not(150<med<2500): continue
      plt.plot(l,y,ls[m],marker='.',color=cols[k%10],lw=1.6,label=f'({m}) {nm[m]} {float(e)*100:g}%: χ²={r["chi2"]:.0f}, r²={r["r2"]:.2f}');k+=1
  plt.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=5,capsize=2,label='Planck TT')
  plt.ylim(0,2600);plt.xlim(1.5,30.5);plt.xlabel('ℓ');plt.ylabel('D_ℓ (μK²)')
  plt.title(f'Realistic web, R={R} Gpc, linear, curves in Planck range (ΛCDM χ²=20.1)')
  plt.legend(loc='lower center',ncol=3,fontsize=8,framealpha=0.85);plt.tight_layout();plt.savefig(f'realweb_R{R}_linP.png',dpi=110);plt.close()
