import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
R=json.load(open('red_inv_t0.5_spectrum.json'));Pl=np.loadtxt('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)]
L=np.arange(2,31);fig,ax=plt.subplots(1,2,figsize=(15,6))
for a,lg in zip(ax,(False,True)):
  a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='tab:blue',ms=4,capsize=2,label='Planck 2018 TT')
  a.plot(L,np.array(R['LCDM']['Dl'])[2:31],'k-',lw=2,label='ΛCDM  χ²=%.1f'%R['LCDM']['chi2'])
  for k,c,s in(('red_inv_t0.5','tab:red','-'),('inv_model_Ag1083','tab:orange','--'),('inv_model_Ag313','tab:green',':')):
    v=R[k];lab={'red_inv_t0.5':'red_inv_t0.5 smooth (= motion term)','inv_model_Ag1083':'motion + Ag·gravity, Ag=1083 (planned)','inv_model_Ag313':'motion + Ag·gravity, Ag=313 (best-fit)'}[k]
    a.plot(L,np.maximum(np.array(v['Dl'])[2:31],1e-3 if lg else -1e9),s,color=c,lw=2,marker='.',label='%s  χ²=%.1f'%(lab,v['chi2']))
  a.set_xlabel('ℓ');a.set_ylabel('D_ℓ [μK²]');a.set_xlim(1.5,30.5)
  if lg:a.set_yscale('log');a.set_title('log y')
  else:a.set_title('linear y')
  a.legend(fontsize=8)
fig.suptitle('Universe red_inv_t0.5 (inversion-prescribed): C=(-0.0392,-0.0501,0.1103,0.025), frac 20637, D/R 3.079')
plt.tight_layout();plt.savefig('red_inv_t0.5_spectrum.png',dpi=130)
