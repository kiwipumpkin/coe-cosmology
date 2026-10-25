import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
S=np.load('mc_sky.npy',mmap_mode='r');sm=np.load('/workspace/redweb2/sm_red.npy');RS=[0.5,0.9,2.0,4.0,8.0,14.0]
l=np.arange(2,31);fig,ax=plt.subplots(3,2,figsize=(16,16));ax=ax.ravel()
for iR,R in enumerate(RS):
  A=[];c=[]
  for s in S:D=Dl(sm+s[iR]);A.append(D[2:31]);c.append(chi(D))
  A=np.array(A);c=np.array(c);o=np.argsort(c);a=ax[iR]
  a.fill_between(l,np.percentile(A,5,0),np.percentile(A,95,0),color='0.8',label=f'5–95% of {len(A)} runs')
  a.plot(l,np.median(A,0),'k-',lw=2,label=f'median run spectrum (median χ²={np.median(c):.0f})')
  for j,k in enumerate(o[:5]):a.plot(l,A[k],'.-',lw=1.2,color=plt.cm.tab10(j),label=f'best #{j+1}: χ²={c[k]:.0f}, r²={r2(np.r_[0,0,A[k]]):.2f}')
  a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='r',ms=4,capsize=2,label='Planck TT')
  a.set_yscale('log');a.set_ylim(1e-1,None);a.set_xlabel('ℓ');a.set_ylabel('D_ℓ (μK²)');a.set_title(f'Monte Carlo, untuned realistic webs, R={R} Gpc (ΛCDM χ²=20.1)')
  a.legend(loc='lower center',ncol=2,fontsize=7,framealpha=.85)
plt.tight_layout();plt.savefig('mc_spectra.png',dpi=95)
