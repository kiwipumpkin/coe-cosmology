import numpy as np
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
t=np.loadtxt('/workspace/joint/lcdm.txt');Lt=np.interp(l,t[:,0],t[:,1])
z=np.load('mc_mean_passing.npz');fig,ax=plt.subplots(3,1,figsize=(11,15),sharex=True)
for a,R in zip(ax,[4.0,8.0,14.0]):
  D=z[f'D{R}'];m=D.mean(0);lo,hi=D.min(0),D.max(0);p16,p84=np.percentile(D,[16,84],0)
  print(f'R={R}: '+'; '.join(f'l={L} min {lo[L]:.0f} max {hi[L]:.0f}' for L in range(2,6)),'| clipped max at l:',[L for L in range(2,31) if hi[L]>2500])
  a.fill_between(l,lo[2:31],hi[2:31],color='lightblue',alpha=.6,label=f'min–max of {len(D)} passing')
  a.fill_between(l,p16[2:31],p84[2:31],color='tab:blue',alpha=.35,label='16–84 percentile')
  a.plot(l,m[2:31],'-',color='tab:blue',lw=2,label=f'mean: χ²={chi(m):.1f}, r²={r2(m):.2f}')
  a.plot(l,Lt,'k-',lw=1.6,label='ΛCDM: χ²=20.1')
  a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',mfc='w',ms=5,capsize=2,label='Planck TT')
  cl=[L for L in range(2,31) if hi[L]>2500]
  a.set_ylim(0,2500);a.set_ylabel('D_ℓ [μK²]');a.set_title(f'R = {R:g} Gpc'+(f'  (max clipped at ℓ={",".join(map(str,cl))})' if cl else ''))
  a.legend(loc='lower center',ncol=3,fontsize=9)
ax[-1].set_xlabel('ℓ');plt.tight_layout();plt.savefig('mc_domain_passing_lin.png',dpi=120)
