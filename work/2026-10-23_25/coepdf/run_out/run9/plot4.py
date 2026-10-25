import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
Pl=np.loadtxt('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)]
l=np.arange(2,31)
def r2(D):return float(np.corrcoef(np.array(D)[2:31],np.interp(l,Pl[:,0],Pl[:,1]))[0,1]**2)
def fig(curves,fn,title):
  plt.figure(figsize=(11,7));cm=plt.cm.viridis(np.linspace(0,1,len(curves)))
  for (lab,D,c2),col in zip(curves,cm):plt.plot(l,np.array(D)[2:31],color='k' if lab.startswith('smooth') else col,lw=2.5 if lab.startswith('smooth') else 1.2,label=f'{lab}: χ²={c2:.4g}, r²={r2(D):.2f}')
  plt.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='r',ms=4,label='Planck TT');plt.yscale('log');plt.xlabel('ℓ');plt.ylabel('D_ℓ (μK²)');plt.title(title)
  plt.legend(fontsize=6.5,ncol=2,loc='upper left',bbox_to_anchor=(1.01,1));plt.tight_layout();plt.savefig(fn,dpi=120);plt.close()
