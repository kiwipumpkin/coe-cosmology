# Fig 25: per column (A | B): overlay of data + CoE (before) + after n=1-3; then residuals after +n=1, +n=1-2, n=1-3.
import sys;sys.path.insert(0,'..')
exec(open('recon_n3.py').read().split("cols=['purple'")[0])
from matplotlib.lines import Line2D
cols=['purple','tab:blue','tab:red','tab:green']
plt.rcParams.update({'font.size':14})
fig=plt.figure(figsize=(16,12.7),facecolor='white')
gs=fig.add_gridspec(4,2,height_ratios=[1.5,1,1,1],hspace=0.45,wspace=0.16,left=.07,right=.99,top=.97,bottom=.105)
box=dict(fc='white',ec='none',alpha=.8)
def bands(x):x.axvspan(15,80,color='gold',alpha=.18,zorder=0);x.axvspan(260,335,color='lightsteelblue',alpha=.35,zorder=0)
for k,(nm,Y,sY,mb,mg) in enumerate((('A',A0,sA,AFl,AFg),('B',B0,sB,BLl,BLg))):
  ms=[(mb.copy(),mg.copy())]
  for n,p in steps2:ms.append((ms[-1][0]+T(n,p,th,k),ms[-1][1]+T(n,p,gt,k)))
  st=lambda m:(np.sum((Y-m)**2/sY**2),np.corrcoef(Y,m)[0,1]**2)
  a=fig.add_subplot(gs[0,k]);bands(a)
  a.errorbar(l,Y,yerr=sY,fmt='o',color='k',ms=3,alpha=.5,lw=.8)
  a.plot(g,ms[0][1],color='purple',lw=2.5);a.plot(g,ms[3][1],color='tab:green',lw=2.5)
  (c0,r0),(c3,r3)=st(ms[0][0]),st(ms[3][0])
  a.text(.02,.05,f'before (CoE): χ²={c0:.2f}, r²={r0:.3f}\nafter n=1–3: χ²={c3:.2f}, r²={r3:.3f}',transform=a.transAxes,fontsize=12.5,bbox=box)
  a.set_title(f'{nm}: CoE model before and after n=1–3',fontweight='bold',fontsize=15)
  a.set_ylabel('A/A$_{all}$' if nm=='A' else 'B/B$_{all}$');a.set_xlim(0,360);a.set_xticks([0,90,180,270,360])
  b0=None
  for j,i in enumerate((1,2,3)):
    b=fig.add_subplot(gs[j+1,k],sharey=b0);b0=b0 or b;bands(b)
    r=Y-ms[i][0];rms=np.sqrt(np.mean(r**2))
    b.errorbar(l,r,yerr=sY,fmt='o',color='k',ms=3,alpha=.5,lw=.8);b.axhline(0,color='grey',lw=1)
    if i<3:
      n2,p2=steps2[i];b.plot(g,p2[2]*np.cos(n2*gt)+p2[3]*np.sin(n2*gt),color=cols[i+1],lw=2.5)
      t=f'RMS={rms:.3f}; next: n={n2}, C={np.hypot(p2[2],p2[3]):.3f}'
    else:t=f'RMS={rms:.3f}; final residual'
    b.text(.02,.06,t,transform=b.transAxes,fontsize=12,bbox=box)
    b.set_title(f'{nm} residual after '+('+n=1' if i==1 else '+n=1–2' if i==2 else 'n=1–3'),fontweight='bold',fontsize=14)
    b.set_ylabel(f'{nm} − model');b.set_xlim(0,360);b.set_xticks([0,90,180,270,360])
    if i==3:b.set_xlabel('ℓ [deg]')
h=[Line2D([],[],color='k',marker='o',ls='',alpha=.5,label='reconstruction, 313 clusters (±σ Migkas Fig. 6)'),Line2D([],[],color='purple',lw=3,label='CoE model, before (A: AF_smooth, B: sheet BL)'),Line2D([],[],color='tab:green',lw=3,label='top: after n=1–3; residuals: next n=3'),Line2D([],[],color='tab:red',lw=3,label='next component n=2'),plt.Rectangle((0,0),1,1,color='gold',alpha=.3,label='ℓ 15–80°'),plt.Rectangle((0,0),1,1,color='lightsteelblue',alpha=.5,label='ℓ 260–335°')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=12.5)
fig.savefig('recon_n3_final.png',dpi=110,facecolor='white')
