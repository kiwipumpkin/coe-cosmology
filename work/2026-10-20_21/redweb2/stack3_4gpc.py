import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
Pl=np.loadtxt(B+'/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)];l=np.arange(2,31)
t=json.load(open('thin4_results.json'));fo=json.load(open('force4_results.json'));dw=json.load(open('densweb_results.json'))
L=t['thin']['0.01'];b=min(L,key=lambda d:d['chi2'])
cur=[('Original web, 1% of bodies kept (best of 5)',b['Dl'],'tab:blue',b['chi2']),
     ('Original web, all bodies at 2% strength',fo['orig']['0.02']['Dl'],'tab:orange',fo['orig']['0.02']['chi2']),
     ('Density-following web (seed 2), 1% strength',dw['s2_R4.0_0.01']['Dl'],'tab:green',dw['s2_R4.0_0.01']['chi2'])]
fig,ax=plt.subplots(figsize=(11,8))
ax.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=5,label='Planck',zorder=10)
for lab,D,c,chi in cur:
  D=np.array(D)[2:31];r2=np.corrcoef(D,Pl[:,1])[0,1]**2
  ax.plot(l,D,color=c,lw=2.2,marker='o',ms=4,label=f'{lab}: χ²={chi:.1f}, r²={r2:.2f}');print(lab,round(chi,1),round(r2,3))
ax.plot(l,np.array(t['smooth']['Dl'])[2:31],'r--',lw=1.8,label=f"Smooth red universe, no web: χ²={t['smooth']['chi2']:.1f}")
ax.set_ylim(0,3000);ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)')
ax.set_title('Best cases, web read at 4 Gpc, vs Planck ℓ=2–30')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.1),ncol=1,fontsize=11,frameon=False)
fig.tight_layout();fig.savefig('stack3_4gpc.png',dpi=110,bbox_inches='tight')
