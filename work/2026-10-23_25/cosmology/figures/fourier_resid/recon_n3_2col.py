# Two-column A | B version: top before/after all harmonics, then harmonic steps stacked.
import sys;sys.path.insert(0,'..')
exec(open('recon_n3.py').read().split("cols=['purple'")[0])
from matplotlib.lines import Line2D
cols=['purple','tab:blue','tab:red','tab:green']
plt.rcParams.update({'font.size':13})
order=[(0,'BEFORE: CoE model, no harmonics'),(3,'AFTER: all harmonics n=1–3 removed'),(1,'step: + shared n=1'),(2,'step: + shared n=1–2 (+n=3 gives AFTER)')]
hr=[]
for _ in order:hr+=[1.25,0.75]
fig,axs=plt.subplots(8,2,figsize=(16,17),facecolor='white',sharex=True,gridspec_kw=dict(height_ratios=hr))
for k,(nm,Y,sY,mb,mg) in enumerate((('A',A0,sA,AFl,AFg),('B',B0,sB,BLl,BLg))):
  ms=[(mb.copy(),mg.copy())]
  for n,p in steps2:ms.append((ms[-1][0]+T(n,p,th,k),ms[-1][1]+T(n,p,gt,k)))
  for row,(i,tt) in enumerate(order):
    m,mgg=ms[i];a=axs[2*row,k];b=axs[2*row+1,k]
    for x in (a,b):x.axvspan(15,80,color='gold',alpha=.18,zorder=0);x.axvspan(260,335,color='lightsteelblue',alpha=.35,zorder=0)
    a.errorbar(l,Y,yerr=sY,fmt='o',color='k',ms=3,alpha=.5,lw=.8)
    if i>0 and row>=2:a.plot(g,ms[i-1][1],color='grey',ls='--',lw=1.6,alpha=.6)
    if row==1:a.plot(g,ms[0][1],color='purple',ls='--',lw=1.6,alpha=.6)
    a.plot(g,mgg,color=cols[i],lw=2.5)
    r=Y-m;c2=np.sum(r**2/sY**2);r2=np.corrcoef(Y,m)[0,1]**2;rms=np.sqrt(np.mean(r**2))
    a.text(.02,.04,f'χ²={c2:.2f}   r²={r2:.3f}',transform=a.transAxes,fontsize=12,bbox=dict(fc='white',ec='none',alpha=.8))
    a.set_title(f'{nm}: {tt}',fontweight='bold',fontsize=13)
    b.errorbar(l,r,yerr=sY,fmt='o',color='k',ms=3,alpha=.5,lw=.8);b.axhline(0,color='grey',lw=1)
    if i<3 and row!=1:
      n2,p2=steps2[i];b.plot(g,p2[2]*np.cos(n2*gt)+p2[3]*np.sin(n2*gt),color=cols[i+1],lw=2.5)
      b.text(.02,.05,f'RMS={rms:.3f}; next: n={n2}, C={np.hypot(p2[2],p2[3]):.3f}',transform=b.transAxes,fontsize=11,bbox=dict(fc='white',ec='none',alpha=.8))
    else:b.text(.02,.05,f'RMS={rms:.3f}; final residual',transform=b.transAxes,fontsize=11,bbox=dict(fc='white',ec='none',alpha=.8))
    a.set_ylabel('A/A$_{all}$' if nm=='A' else 'B/B$_{all}$');b.set_ylabel('resid.')
    a.set_xlim(0,360);a.set_xticks([0,90,180,270,360])
  for j in range(8):
    if j>1:axs[j,k].sharey(axs[j%2,k])
  axs[7,k].set_xlabel('ℓ [deg]')
h=[Line2D([],[],color='k',marker='o',ls='',alpha=.5,label='reconstruction, 313 clusters (±σ Migkas Fig. 6)'),Line2D([],[],color='purple',lw=3,label='CoE model (A: AF_smooth, B: sheet BL)')]+[Line2D([],[],color=cols[j],lw=3,label=f'shared n={j} (model / component)') for j in (1,2,3)]+[Line2D([],[],color='grey',ls='--',lw=2,label='previous model'),plt.Rectangle((0,0),1,1,color='gold',alpha=.3,label='ℓ 15–80°'),plt.Rectangle((0,0),1,1,color='lightsteelblue',alpha=.5,label='ℓ 260–335°')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=11)
fig.tight_layout(rect=(0,.05,1,1),h_pad=0.4,w_pad=1);fig.savefig('recon_n3_2col.png',dpi=110,facecolor='white')
