import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('plot4.py').read());T=json.load(open('thin4_results.json'))
plt.figure(figsize=(11,7));l=np.arange(2,31)
for f,c,m in zip(['0.01','0.0125','0.015','0.0175','0.02'],['tab:blue','tab:green','tab:purple','tab:brown','tab:orange'],'osD^v'):
  L=T['thin'][f];cs=[d['chi2'] for d in L];b=L[int(np.argmin(cs))];print(f,'%.1f %.4g %.2f'%(min(cs),np.mean(cs),r2(b['Dl'])))
  plt.plot(l,np.array(b['Dl'])[2:31],color=c,marker=m,ms=5,lw=1.5,label=f'kept {float(f)*100:g}%: χ²={min(cs):.1f}, r²={r2(b["Dl"]):.2f}')
s=T['smooth'];plt.plot(l,np.array(s['Dl'])[2:31],'r--',lw=2,label=f'smooth red: χ²={s["chi2"]:.1f}, r²={r2(s["Dl"]):.2f}')
plt.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=4,capsize=2,label='Planck TT')
plt.ylim(0,3000);plt.xlabel('ℓ');plt.ylabel('D_ℓ (μK²)');plt.title('Density-shaped thinning, web read at 4 Gpc, best of 5 draws (χ², Pearson r² over ℓ=2–30)');plt.legend();plt.tight_layout();plt.savefig('thin4_linear3.png',dpi=120)
