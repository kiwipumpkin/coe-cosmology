import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
S=np.load('mc_sky.npy');N=len(S);sm=np.load('/workspace/redweb2/sm_red.npy');RS=[0.5,0.9,2.0,4.0,8.0,14.0]
t=np.loadtxt('/workspace/joint/lcdm.txt');Dth=np.zeros(31);Dth[2:31]=np.interp(np.arange(2,31),t[:,0],t[:,1])
rng=np.random.default_rng(7);L=np.arange(31);LC=[];LR=[]
for i in range(1000):
  D=np.zeros(31);D[2:]=Dth[2:]*rng.chisquare(2*L[2:]+1)/(2*L[2:]+1);LC.append(chi(D));LR.append(r2(D))
LC=np.array(LC);LR=np.array(LR)
def summ(c,r):return dict(N=len(c),chi2_median=float(np.median(c)),chi2_p5=float(np.percentile(c,5)),chi2_p95=float(np.percentile(c,95)),frac_le_20p1=float(np.mean(c<=20.1)),frac_le_42p6=float(np.mean(c<=42.6)),r2_median=float(np.median(r)),r2_p5=float(np.percentile(r,5)),r2_p95=float(np.percentile(r,95)))
res=dict(lcdm=summ(LC,LR),realweb={})
fig,ax=plt.subplots(2,3,figsize=(17,9));ax=ax.ravel()
for iR,R in enumerate(RS):
  c=[];r=[]
  for s in S:D=Dl(sm+s[iR]);c.append(chi(D));r.append(r2(D))
  c=np.array(c);r=np.array(r);d=summ(c,r)
  for k in ('20p1','42p6'):a=d['frac_le_'+k];b=res['lcdm']['frac_le_'+k];d['ratio_frac_le_'+k+'_vs_lcdm']=a/b if b>0 else None
  res['realweb'][str(R)]=d;np.save(f'mc_chi2_R{R}.npy',c)
  bins=np.logspace(np.log10(min(c.min(),LC.min())*0.8),np.log10(max(c.max(),LC.max())*1.2),60)
  ax[iR].hist(LC,bins,alpha=.6,label=f'ΛCDM cosmic variance (N=1000), median {np.median(LC):.1f}');ax[iR].hist(c,bins,alpha=.6,label=f'realistic red web (N={N}), median {np.median(c):.1f}')
  for v,ls in ((20.1,':'),(42.6,'--')):ax[iR].axvline(v,color='k',ls=ls)
  ax[iR].set_xscale('log');ax[iR].set_title(f'R={R} Gpc'+('  (nearest last scattering ~14 Gpc)' if R==14 else ''));ax[iR].set_xlabel('χ² vs Planck ℓ=2–30');ax[iR].legend(fontsize=7)
plt.tight_layout();plt.savefig('mc_chi2_hist.png',dpi=105);json.dump(res,open('mc_results.json','w'),indent=1)
print(json.dumps(res,indent=0))
