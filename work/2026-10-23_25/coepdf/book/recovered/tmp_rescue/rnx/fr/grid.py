import numpy as np
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':17})
d=np.genfromtxt('residuals.csv',delimiter=',',names=True)
l,r,s=d['ell_deg'],d['residual'],d['sigma'];w=1/s**2;th=np.radians(l)
D=lambda t,n:np.column_stack([np.ones_like(t)]+[f(k*t) for k in range(1,n+1) for f in(np.cos,np.sin)])
fit=lambda y:np.linalg.lstsq(D(th,6)*np.sqrt(w)[:,None],y*np.sqrt(w),rcond=None)[0]
comp=lambda p,n,t:p[2*n-1]*np.cos(n*t)+p[2*n]*np.sin(n*t)
chi=lambda y:float(np.sum(w*y**2))
g=np.linspace(0,360,721);gt=np.radians(g);cols=['tab:blue','tab:red','tab:green','tab:orange']
fig,axs=plt.subplots(5,2,figsize=(20,26),sharex=True,sharey=True)
pA=fit(r);yA=r.copy();yB=r.copy();rows=[]
for k in range(5):
 if k>0:
  yA=yA-comp(pA,k,th); pB=fit(yB); Ck=np.hypot(pB[2*k-1],pB[2*k]); phk=np.degrees(np.arctan2(pB[2*k],pB[2*k-1]))%360/k
  yB=yB-comp(pB,k,th); rows.append((k,Ck,phk))
 for j,(ax,y,p) in enumerate(((axs[k,0],yA,pA),(axs[k,1],yB,fit(yB)))):
  ax.errorbar(l,y,yerr=s,fmt='o',color='k',ms=4,alpha=.6,capsize=2)
  tot=np.full_like(g,p[0])
  for n in range(1,5):
   if j==0 and n<=k: continue
   c=comp(p,n,gt);tot+=c;ax.plot(g,c,color=cols[n-1],lw=2.5)
  ax.plot(g,tot,'k--',lw=2.2);ax.axhline(0,color='grey',lw=1)
  rem='none' if k==0 else 'n='+','.join(map(str,range(1,k+1)))
  ax.set_title(f"{'A full-fit' if j==0 else 'B iterative'} | removed: {rem} | χ²={chi(y):.2f}",fontsize=17)
  print(j,k,round(chi(y),2))
 axs[k,0].set_ylabel('residual')
for ax in axs[-1]:ax.set_xlabel('Galactic longitude ℓ [deg]')
axs[0,0].set_xlim(0,360)
from matplotlib.lines import Line2D
h=[Line2D([],[],color=c,lw=3,label=f'component n={i+1}') for i,c in enumerate(cols)]+[Line2D([],[],color='k',ls='--',lw=2.5,label='sum of shown terms + c₀'),Line2D([],[],color='k',marker='o',ls='',label='residual ±σ')]
fig.legend(handles=h,loc='lower center',ncol=3,fontsize=19,frameon=False)
fig.suptitle('Axisymmetric A/A_all residuals: removing harmonics one at a time (A: from single n≤6 fit; B: refit each step)',fontsize=20)
fig.tight_layout(rect=(0,0.05,1,0.98));fig.savefig('harmonic_removal_grid_A_vs_B.png',dpi=90)
for x in rows:print('B',x)
