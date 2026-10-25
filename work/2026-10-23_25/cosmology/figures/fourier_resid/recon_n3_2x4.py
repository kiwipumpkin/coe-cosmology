# Term-by-term shared harmonic decomposition of A and B (chain: cum2.py -> recon_n3.py), with residuals beneath.
import sys;sys.path.insert(0,'..')
exec(open('recon_n3.py').read().split("cols=['purple'")[0])
from matplotlib.lines import Line2D
cols=['purple','tab:blue','tab:red','tab:green']
plt.rcParams.update({'font.size':19})
fig,axs=plt.subplots(4,4,figsize=(26,21),facecolor='white',sharex=True,gridspec_kw=dict(height_ratios=[1.25,0.8,1.25,0.8]))
stats={}
for blk,(nm,Y,sY,mb,mg,k) in enumerate((('A',A0,sA,AFl,AFg,0),('B',B0,sB,BLl,BLg,1))):
  m,mgg=mb.copy(),mg.copy()
  for i in range(4):
    if i>0:n,p=steps2[i-1];pm=mgg.copy();m=m+T(n,p,th,k);mgg=mgg+T(n,p,gt,k)
    a=axs[2*blk,i];b=axs[2*blk+1,i]
    for x in (a,b):x.axvspan(15,80,color='gold',alpha=.18,zorder=0);x.axvspan(260,335,color='lightsteelblue',alpha=.35,zorder=0)
    a.errorbar(l,Y,yerr=sY,fmt='o',color='k',ms=4,alpha=.5,lw=1)
    if i>0:a.plot(g,pm,color='grey',ls='--',lw=2,alpha=.6)
    a.plot(g,mgg,color=cols[i],lw=3.5)
    r=Y-m;c2=np.sum(r**2/sY**2);r2=np.corrcoef(Y,m)[0,1]**2;rms=np.sqrt(np.mean(r**2));stats[(nm,i)]=(rms,c2,r2)
    a.text(.03,.04,f'χ²={c2:.2f}   r²={r2:.3f}',transform=a.transAxes,fontsize=19,color='0.15',bbox=dict(fc='white',ec='none',alpha=.8))
    b.errorbar(l,r,yerr=sY,fmt='o',color='k',ms=4,alpha=.5,lw=1);b.axhline(0,color='grey',lw=1.2)
    if i<3:
      n2,p2=steps2[i];b.plot(g,p2[2]*np.cos(n2*gt)+p2[3]*np.sin(n2*gt),color=cols[i+1],lw=3.5)
      b.text(.03,.05,f'RMS={rms:.3f}; next: n={n2}, C={np.hypot(p2[2],p2[3]):.3f}',transform=b.transAxes,fontsize=17,color='0.15',bbox=dict(fc='white',ec='none',alpha=.8))
    else:b.text(.03,.05,f'RMS={rms:.3f}; final residual',transform=b.transAxes,fontsize=17,color='0.15',bbox=dict(fc='white',ec='none',alpha=.8))
    if True:a.set_title('CoE model (before)' if i==0 else ('+ shared n=1' if i==1 else ('+ shared n=1–2' if i==2 else 'after all: n=1–3')),fontweight='bold',fontsize=21)
    if i==0:a.text(-.32,.5,nm,transform=a.transAxes,fontsize=34,fontweight='bold',va='center')
    a.set_xlim(0,360);a.set_xticks([0,90,180,270,360])
    if i>0:a.sharey(axs[2*blk,0]);b.sharey(axs[2*blk+1,0]);a.tick_params(labelleft=False);b.tick_params(labelleft=False)
  axs[2*blk,0].set_ylabel('A/A$_{all}$' if nm=='A' else 'B/B$_{all}$',fontsize=22);axs[2*blk+1,0].set_ylabel(f'{nm} residual',fontsize=22)
for b in axs[3]:b.set_xlabel('ℓ [deg]',fontsize=21)
h=[Line2D([],[],color='k',marker='o',ls='',alpha=.5,label='reconstruction, 313 clusters (±σ Migkas Fig. 6)'),Line2D([],[],color='purple',lw=3,label='CoE model (A: AF_smooth, B: sheet BL)')]+[Line2D([],[],color=cols[j],lw=3,label=f'shared n={j} (model / component)') for j in (1,2,3)]+[Line2D([],[],color='grey',ls='--',lw=2,label='previous model'),plt.Rectangle((0,0),1,1,color='gold',alpha=.3,label='ℓ 15–80°'),plt.Rectangle((0,0),1,1,color='lightsteelblue',alpha=.5,label='ℓ 260–335°')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=19)
fig.tight_layout(rect=(0,.075,1,1),h_pad=0.6,w_pad=0.4);fig.savefig('recon_n3_2x4.png',dpi=90,facecolor='white')
for kk,v in stats.items():print(kk,'RMS=%.3f chi2=%.2f r2=%.3f'%v)
for n,p in steps2:print('n',n,'C=%.4f'%np.hypot(p[2],p[3]),'phase=%.1f'%(np.degrees(np.arctan2(p[3],p[2]))%360),'offA=%.4f offB=%.4f'%(p[0],p[1]))
