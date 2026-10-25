# Fig 25: 2 columns (A | B) x 4 model rows (CoE, +n=1, +n=1-2, after all n=1-3); residual strip under row 1 only.
import sys;sys.path.insert(0,'..')
exec(open('recon_n3.py').read().split("cols=['purple'")[0])
from matplotlib.lines import Line2D
cols=['purple','tab:blue','tab:red','tab:green']
titles=['CoE model (before)','+ shared n=1','+ shared n=1–2','after all: n=1–3']
plt.rcParams.update({'font.size':14})
fig=plt.figure(figsize=(16,16.5),facecolor='white')
gs=fig.add_gridspec(5,2,height_ratios=[1,0.45,1,1,1],hspace=0.42,wspace=0.16,left=.07,right=.99,top=.975,bottom=.1)
for k,(nm,Y,sY,mb,mg) in enumerate((('A',A0,sA,AFl,AFg),('B',B0,sB,BLl,BLg))):
  ms=[(mb.copy(),mg.copy())]
  for n,p in steps2:ms.append((ms[-1][0]+T(n,p,th,k),ms[-1][1]+T(n,p,gt,k)))
  a0=None
  for i in range(4):
    a=fig.add_subplot(gs[0 if i==0 else i+1,k],sharey=a0);a0=a0 or a
    m,mgg=ms[i]
    a.axvspan(15,80,color='gold',alpha=.18,zorder=0);a.axvspan(260,335,color='lightsteelblue',alpha=.35,zorder=0)
    a.errorbar(l,Y,yerr=sY,fmt='o',color='k',ms=3,alpha=.5,lw=.8)
    if i>0:a.plot(g,ms[i-1][1],color='grey',ls='--',lw=1.6,alpha=.6)
    a.plot(g,mgg,color=cols[i],lw=2.5)
    r=Y-m;c2=np.sum(r**2/sY**2);r2=np.corrcoef(Y,m)[0,1]**2
    a.text(.02,.05,f'χ²={c2:.2f}   r²={r2:.3f}',transform=a.transAxes,fontsize=13,bbox=dict(fc='white',ec='none',alpha=.8))
    a.set_title(f'{nm}: {titles[i]}',fontweight='bold',fontsize=15)
    a.set_ylabel('A/A$_{all}$' if nm=='A' else 'B/B$_{all}$');a.set_xlim(0,360);a.set_xticks([0,90,180,270,360])
    if i==3:a.set_xlabel('ℓ [deg]')
    if i==0:
      b=fig.add_subplot(gs[1,k],sharex=a);rms=np.sqrt(np.mean(r**2))
      b.axvspan(15,80,color='gold',alpha=.18,zorder=0);b.axvspan(260,335,color='lightsteelblue',alpha=.35,zorder=0)
      b.errorbar(l,r,yerr=sY,fmt='o',color='k',ms=3,alpha=.5,lw=.8);b.axhline(0,color='grey',lw=1)
      n2,p2=steps2[0];b.plot(g,p2[2]*np.cos(n2*gt)+p2[3]*np.sin(n2*gt),color=cols[1],lw=2.5)
      b.text(.02,.06,f'RMS={rms:.3f}; next: n=1, C={np.hypot(p2[2],p2[3]):.3f}',transform=b.transAxes,fontsize=12,bbox=dict(fc='white',ec='none',alpha=.8))
      b.set_ylabel(f'{nm} residual')
h=[Line2D([],[],color='k',marker='o',ls='',alpha=.5,label='reconstruction, 313 clusters (±σ Migkas Fig. 6)'),Line2D([],[],color='purple',lw=3,label='CoE model (A: AF_smooth, B: sheet BL)')]+[Line2D([],[],color=cols[j],lw=3,label=f'shared n={j} (model / component)') for j in (1,2,3)]+[Line2D([],[],color='grey',ls='--',lw=2,label='previous model'),plt.Rectangle((0,0),1,1,color='gold',alpha=.3,label='ℓ 15–80°'),plt.Rectangle((0,0),1,1,color='lightsteelblue',alpha=.5,label='ℓ 260–335°')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=13)
fig.savefig('recon_n3_2c4r.png',dpi=110,facecolor='white')
