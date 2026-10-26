# PROVENANCE COPY of /workspace/cosmology/figures/fourier_resid/recon_n3_2x2.py (original generating script, rev6 figure; run in its original directory).
# Output figure copied here as L25_recon_n3.png. Not re-run for rev7: see Appendix B.
# Fig 25 2x2: top = data + CoE before + after n=1-3; bottom = residuals after +n=1, +n=1-2, n=1-3 overlaid.
import sys;sys.path.insert(0,'/workspace/cosmology/figures')
exec(open('recon_n3.py').read().split("cols=['purple'")[0])
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size':14})
RC=['tab:blue','tab:orange','magenta'];lab=['after +n=1','after +n=1–2','after n=1–3']
fig=plt.figure(figsize=(16,11),facecolor='white')
gs=fig.add_gridspec(2,2,height_ratios=[1,1],hspace=0.32,wspace=0.16,left=.07,right=.99,top=.95,bottom=.24)
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
  a.set_title(f'{nm}: data, CoE model before and after n=1–3',fontweight='bold',fontsize=15)
  a.set_ylabel('A/A$_{all}$' if nm=='A' else 'B/B$_{all}$');a.set_xlim(0,360);a.set_xticks([0,90,180,270,360])
  b=fig.add_subplot(gs[1,k]);bands(b);b.axhline(0,color='grey',lw=1)
  hh=[]
  for j,i in enumerate((1,2,3)):
    r=Y-ms[i][0];rms=np.sqrt(np.mean(r**2))
    if i==3:b.errorbar(l,r,yerr=sY,fmt='none',ecolor='0.6',lw=.7,alpha=.6,zorder=1)
    hh+=b.plot(l,r,':o',color=RC[j],ms=3.5,lw=0.9,alpha=.8,label=f'{lab[j]}: RMS={rms:.3f}')
    hh+=b.plot(g,T(*steps2[j],gt,k),'-',color=RC[j],lw=2.6,label=f'fitted n={i} mode')
  b.legend(handles=hh,loc='upper center',bbox_to_anchor=(.5,-.2),ncol=2,fontsize=11,frameon=False)
  b.set_title(f'{nm} residuals (data − model), overlaid',fontweight='bold',fontsize=15)
  b.set_ylabel(f'{nm} − model');b.set_xlim(0,360);b.set_xticks([0,90,180,270,360]);b.set_xlabel('ℓ [deg]')
h=[Line2D([],[],color='k',marker='o',ls='',alpha=.5,label='reconstruction, 313 clusters (±σ Migkas Fig. 6; grey bars on n=1–3 residual)'),Line2D([],[],color='purple',lw=3,label='CoE model, before (A: AF_smooth, B: sheet BL)'),Line2D([],[],color='tab:green',lw=3,label='after shared n=1–3 removed')]+[Line2D([],[],color=RC[j],marker='o',lw=1.4,label='residual '+lab[j]) for j in range(3)]+[plt.Rectangle((0,0),1,1,color='gold',alpha=.3,label='ℓ 15–80°'),plt.Rectangle((0,0),1,1,color='lightsteelblue',alpha=.5,label='ℓ 260–335°')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=12.5)
fig.savefig('/workspace/coepdf/fig7/rerun/recon_n3_2x2.png',dpi=110,facecolor='white')
