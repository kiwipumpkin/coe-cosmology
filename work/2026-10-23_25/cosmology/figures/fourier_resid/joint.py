import numpy as np,csv
from scipy.stats import pearsonr
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':17})
from matplotlib.lines import Line2D
rows=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
l=np.array([float(r['l_deg']) for r in rows]);th=np.radians(l)
d=np.genfromtxt('residuals.csv',delimiter=',',names=True);rA,sA=d['residual'],d['sigma']
Bv=np.array([float(r['B_over_Ball']) for r in rows]);sB=np.array([(float(r['B_err_low'])+float(r['B_err_high']))/2 for r in rows]);wB=1/sB**2;rB=Bv-np.sum(wB*Bv)/np.sum(wB)
H=lambda t,N:np.column_stack([f(k*t) for k in range(1,N+1) for f in(np.cos,np.sin)]) if N else np.zeros((len(t),0))
def solve(N,k):
 n=len(l);Z=np.zeros(n);o=np.ones(n);h=H(th,N)
 X=np.vstack([np.column_stack([o,Z,h]),np.column_stack([Z,o,k*h])]);y=np.r_[rA,rB];s=np.r_[sA,sB]
 p=np.linalg.lstsq(X/s[:,None],y/s,rcond=None)[0];m=X@p
 cA=np.sum((rA-m[:n])**2/sA**2);cB=np.sum((rB-m[n:])**2/sB**2);return p,cA,cB,m[:n],m[n:]
res={};print('N variant k chiA chiB comb r2A r2B Cn/peak')
for N in range(0,5):
 for v,ks in(('+',[1.0]),('-',[-1.0]),('k',np.linspace(-3,3,601))):
  best=min(((solve(N,k),k) for k in ks),key=lambda z:z[0][1]+z[0][2])
  (p,cA,cB,mA,mB),k=best;res[(N,v)]=(p,cA,cB,k)
  r2A=pearsonr(rA,mA)[0]**2 if N else np.nan;r2B=pearsonr(rB,mB)[0]**2 if N else np.nan
  cn=' '.join(f'n{i}:{np.hypot(p[2*i],p[2*i+1]):.3f}@{np.degrees(np.arctan2(p[2*i+1],p[2*i]))%360/i:.0f}' for i in range(1,N+1))
  print(N,v,round(k,2),round(cA,2),round(cB,2),round(cA+cB,2),round(r2A,3),round(r2B,3),cn)
# simple average of separate fits (n<=4 terms from each full n<=6 fit)
def sep(r,s):
 X=np.column_stack([np.ones_like(th),H(th,6)]);return np.linalg.lstsq(X/s[:,None],r/s,rcond=None)[0]
pa,pb=sep(rA,sA),sep(rB,sB)
for N in range(1,5):
 avg=(pa[1:2*N+1]+pb[1:2*N+1])/2;m=H(th,N)@avg
 cA=np.sum((rA-pa[0]-m)**2/sA**2);cB=np.sum((rB-pb[0]-m)**2/sB**2);print('avg',N,round(cA,2),round(cB,2),round(cA+cB,2))
np.save('joint_res.npy',np.array([res],dtype=object))
g=np.linspace(0,360,721);gt=np.radians(g);cols=['tab:blue','tab:red','tab:green','tab:orange']
fig,axs=plt.subplots(5,2,figsize=(20,26),sharex=True,sharey=True)
for N in range(5):
 p,cA,cB,k=res[(N,'-')]
 for j,(r,s,c0,sg,c2,nm) in enumerate(((rA,sA,p[0],1,cA,'A'),(rB,sB,p[1],-1,cB,'B'))):
  ax=axs[N,j];ax.errorbar(l,r,yerr=s,fmt='o',color='k',ms=4,alpha=.6,capsize=2);tot=np.full_like(g,c0)
  for i in range(1,N+1):cc=sg*(p[2*i]*np.cos(i*gt)+p[2*i+1]*np.sin(i*gt));tot+=cc;ax.plot(g,cc,color=cols[i-1],lw=2.5)
  ax.plot(g,tot,'k--',lw=2.2);ax.axhline(0,color='grey',lw=1)
  ax.set_title(f"{nm} | shared harmonics N={N}{'' if j==0 else ' (×−1)'} | χ²={c2:.2f}")
 axs[N,0].set_ylabel('residual')
axs[0,0].text(.5,1.22,'A/A_all − axisymmetric model (AF_smooth)',transform=axs[0,0].transAxes,ha='center',fontsize=19,weight='bold')
axs[0,1].text(.5,1.22,'B/B_all − weighted mean (no B model exists)',transform=axs[0,1].transAxes,ha='center',fontsize=19,weight='bold')
for ax in axs[-1]:ax.set_xlabel('Galactic longitude ℓ [deg]')
axs[0,0].set_xlim(0,360)
h=[Line2D([],[],color=c,lw=3,label=f'shared component n={i+1}') for i,c in enumerate(cols)]+[Line2D([],[],color='k',ls='--',lw=2.5,label='joint model: c₀ + Σ shared (B uses −1×)'),Line2D([],[],color='k',marker='o',ls='',label='residual data ±σ')]
fig.legend(handles=h,loc='lower center',ncol=3,fUcontsize=19,frameon=False) if False else fig.legend(handles=h,loc='lower center',ncol=3,fontsize=19,frameon=False)
fig.suptitle('Joint fit: one shared harmonic set minimising χ²_A + χ²_B (1/σ² weights), B = −shared\nrows: shown residuals with the N-term joint model overlaid (not subtracted)',fontsize=19)
fig.tight_layout(rect=(0,0.05,1,0.955));fig.savefig('joint_shared_harmonics_signflip.png',dpi=90)
