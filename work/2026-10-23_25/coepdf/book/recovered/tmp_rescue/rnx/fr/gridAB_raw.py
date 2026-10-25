import numpy as np,csv
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':17})
from matplotlib.lines import Line2D
rows=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
l=np.array([float(r['l_deg']) for r in rows]);th=np.radians(l)
Bv=np.array([float(r['B_over_Ball']) for r in rows]);sB=np.array([(float(r['B_err_low'])+float(r['B_err_high']))/2 for r in rows])
Av=np.array([float(r['A_over_Aall']) for r in rows]);sA=np.array([(float(r['A_err_low'])+float(r['A_err_high']))/2 for r in rows]);wA=1/sA**2;Amean=np.sum(wA*Av)/np.sum(wA);rA=Av-Amean
wB=1/sB**2;Bmean=np.sum(wB*Bv)/np.sum(wB);rB=Bv-Bmean
D=lambda t,n:np.column_stack([np.ones_like(t)]+[f(k*t) for k in range(1,n+1) for f in(np.cos,np.sin)])
def fit(y,s):w=1/s**2;return np.linalg.lstsq(D(th,6)*np.sqrt(w)[:,None],y*np.sqrt(w),rcond=None)[0]
comp=lambda p,n,t:p[2*n-1]*np.cos(n*t)+p[2*n]*np.sin(n*t)
g=np.linspace(0,360,721);gt=np.radians(g);cols=['tab:blue','tab:red','tab:green','tab:orange']
fig,axs=plt.subplots(5,2,figsize=(20,26),sharex=True,sharey='col')
out={}
for j,(r,s,name) in enumerate(((rA,sA,f'A/A_all − weighted mean ({Amean:.4f})'),(rB,sB,f'B/B_all − weighted mean ({Bmean:.4f})'))):
 p=fit(r,s);y=r.copy();chis=[];cn=[]
 for k in range(5):
  if k>0:y=y-comp(p,k,th)
  c2=float(np.sum(y**2/s**2));chis.append(c2);ax=axs[k,j]
  ax.errorbar(l,y,yerr=s,fmt='o',color='k',ms=4,alpha=.6,capsize=2);tot=np.full_like(g,p[0])
  for n in range(k+1,5):cc=comp(p,n,gt);tot+=cc;ax.plot(g,cc,color=cols[n-1],lw=2.5)
  ax.plot(g,tot,'k--',lw=2.2);ax.axhline(0,color='grey',lw=1)
  rem='none' if k==0 else 'n=1' if k==1 else f'n=1–{k}'
  ax.set_title(f"{'A' if j==0 else 'B'} | removed: {rem} | χ²={c2:.2f}")
 for n in range(1,7):cn.append((n,np.hypot(p[2*n-1],p[2*n]),np.degrees(np.arctan2(p[2*n],p[2*n-1]))%360/n))
 out[name]=(chis,cn);axs[0,j].text(0.5,1.22,name,transform=axs[0,j].transAxes,ha='center',fontsize=19,weight='bold')
for k in range(5):axs[k,0].set_ylabel('residual')
for ax in axs[-1]:ax.set_xlabel('Galactic longitude ℓ [deg]')
axs[0,0].set_xlim(0,360)
h=[Line2D([],[],color=c,lw=3,label=f'component n={i+1}') for i,c in enumerate(cols)]+[Line2D([],[],color='k',ls='--',lw=2.5,label='sum of remaining n≤4 + c₀'),Line2D([],[],color='k',marker='o',ls='',label='residual ±σ (Migkas Fig. 6, 71 bins)')]
fig.legend(handles=h,loc='lower center',ncol=3,fontsize=19,frameon=False)
fig.suptitle('Removing harmonics one at a time (single n≤6 weighted Fourier fit, 360° period)\nBoth columns: raw Migkas data minus its own 1/σ²-weighted mean (no model subtracted)',fontsize=19)
fig.tight_layout(rect=(0,0.05,1,0.955));fig.savefig('harmonic_removal_grid_A_and_B_raw.png',dpi=90)
for k,v in out.items():print(k,[round(x,2) for x in v[0]]);[print('  n=%d C=%.4f peak=%.1f'%c) for c in v[1]]
# reference: AN mass-proxy strip
