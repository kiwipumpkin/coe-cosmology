import numpy as np,csv
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':17})
from matplotlib.lines import Line2D
rows=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
l=np.array([float(r['l_deg']) for r in rows]);th=np.radians(l)
def get(k):
 v=np.array([float(r[f'{k}_over_{k}all']) for r in rows]);s=np.array([(float(r[f'{k}_err_low'])+float(r[f'{k}_err_high']))/2 for r in rows]);w=1/s**2
 return v-np.sum(w*v)/np.sum(w),s
rA,sA=get('A');rB,sB=get('B')
D=lambda t,n:np.column_stack([np.ones_like(t)]+[f(k*t) for k in range(1,n+1) for f in(np.cos,np.sin)])
def fit(y,s):w=1/s**2;return np.linalg.lstsq(D(th,6)*np.sqrt(w)[:,None],y*np.sqrt(w),rcond=None)[0]
comp=lambda p,n,t:p[2*n-1]*np.cos(n*t)+p[2*n]*np.sin(n*t)
g=np.linspace(0,360,721);gt=np.radians(g);cols=['tab:blue','tab:red','tab:green','tab:orange']
for drv,p in(('A',fit(rA,sA)),('B',fit(rB,sB))):
 fig,axs=plt.subplots(5,2,figsize=(20,26),sharex=True,sharey=True);res={}
 for j,(r,s,nm) in enumerate(((rA,sA,'A'),(rB,sB,'B'))):
  y=r.copy();ch=[]
  for k in range(5):
   if k>0:y=y-comp(p,k,th)
   c2=float(np.sum(y**2/s**2));ch.append(c2);ax=axs[k,j]
   ax.errorbar(l,y,yerr=s,fmt='o',color='k',ms=4,alpha=.6,capsize=2);tot=np.full_like(g,p[0])
   for n in range(k+1,5):cc=comp(p,n,gt);tot+=cc;ax.plot(g,cc,color=cols[n-1],lw=2.5)
   ax.plot(g,tot,'k--',lw=2.2);ax.axhline(0,color='grey',lw=1)
   rem='none' if k==0 else f'{drv}-fit n=1' if k==1 else f'{drv}-fit n=1–{k}'
   ax.set_title(f'{nm} | removed: {rem} | χ²={c2:.2f}')
  res[nm]=ch;axs[0,j].text(.5,1.22,f'{nm}/{nm}_all − weighted mean',transform=axs[0,j].transAxes,ha='center',fontsize=19,weight='bold')
 for k in range(5):axs[k,0].set_ylabel('residual')
 for ax in axs[-1]:ax.set_xlabel('Galactic longitude ℓ [deg]')
 axs[0,0].set_xlim(0,360)
 h=[Line2D([],[],color=c,lw=3,label=f'{drv}-fit component n={i+1}') for i,c in enumerate(cols)]+[Line2D([],[],color='k',ls='--',lw=2.5,label=f'{drv}-fit: remaining n≤4 + c₀'),Line2D([],[],color='k',marker='o',ls='',label='residual ±σ (Migkas Fig. 6)')]
 fig.legend(handles=h,loc='lower center',ncol=3,fontsize=19,frameon=False)
 fig.suptitle(f'{drv}-driven: harmonics fitted to raw {drv} (n≤6, 1/σ² weights) are subtracted from BOTH A and B\nNothing fitted to the other observable; χ² = residual about zero',fontsize=19)
 fig.tight_layout(rect=(0,0.05,1,0.955));fig.savefig(f'harmonic_removal_grid_shared_{drv}driven.png',dpi=90)
 print(drv,{k:[round(x,2) for x in v] for k,v in res.items()})
