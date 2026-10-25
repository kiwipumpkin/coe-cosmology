import numpy as np
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
t=np.loadtxt('/workspace/joint/lcdm.txt');Lt=np.interp(l,t[:,0],t[:,1])
z=np.load('mc_mean_passing.npz');fig,ax=plt.subplots(2,1,figsize=(11,12),sharex=True)
X=ax[0];allD=[]
for R,c in zip([4.0,8.0,14.0],['tab:blue','tab:orange','tab:green']):
  D=z[f'D{R}'];allD.append(D);m=D.mean(0);p16,p84=np.percentile(D,[16,84],0)
  X.fill_between(l,D.min(0)[2:31],D.max(0)[2:31],color=c,alpha=.2,lw=0,label=f'R={R:g} min–max (N={len(D)})')
  X.fill_between(l,p16[2:31],p84[2:31],color=c,alpha=.3,lw=0,label=f'R={R:g} 16–84%')
  X.plot(l,m[2:31],color=c,lw=2,label=f'R={R:g} mean: χ²={chi(m):.1f}, r²={r2(m):.2f}')
for A in ax:
  A.plot(l,Lt,'k-',lw=1.6,label='ΛCDM: χ²=20.1')
  A.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',mfc='w',ms=5,capsize=2,label='Planck TT')
  A.set_ylim(0,2500);A.set_ylabel('D_ℓ [μK²]')
X.set_title('Passing-universe domains at R = 4, 8, 14 Gpc overlaid (overlaps darker; ℓ=3 and R=8 ℓ 8–12 max clipped)',fontsize=10)
X.legend(loc='lower center',ncol=3,fontsize=8)
A=np.vstack(allD);yb=np.linspace(0,2500,251);yc=(yb[1:]+yb[:-1])/2
F=np.array([[np.mean((A[:,L]<=y)&False) for y in yc] for L in l])  # placeholder
F=np.array([np.histogram(A[:,L],bins=yb)[0] for L in l],float)
# fraction of all spectra whose min-max... use per-l density of all passing spectra
F/=len(A)
# count of R min-max bands covering
C=np.array([[sum((D[:,L].min()<=y<=D[:,L].max()) for D in allD) for y in yc] for L in l])
im=ax[1].pcolormesh(np.r_[l-.5,30.5],yb,F.T,cmap='Blues',shading='flat',vmin=0)
ax[1].contour(l,yc,C.T,levels=[1.5,2.5],colors=['gray','k'],linewidths=[1,1.5])
plt.colorbar(im,ax=ax[1],pad=.01,label='fraction of all passing spectra (all R) per 10 μK² bin')
ax[1].plot([],[],color='gray',label='contour: ≥2 R min–max bands overlap');ax[1].plot([],[],color='k',lw=1.5,label='contour: all 3 R bands overlap')
ax[1].set_title(f'Overlap density: all {len(A)} passing spectra pooled across R',fontsize=10)
ax[1].legend(loc='lower center',ncol=2,fontsize=8);ax[1].set_xlabel('ℓ')
plt.tight_layout();plt.savefig('mc_domain_overlay_lin.png',dpi=120)
