import numpy as np
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':16})
d=np.genfromtxt('residuals.csv',delimiter=',',names=True)
l,r,s=d['ell_deg'],d['residual'],d['sigma'];w=1/s**2;th=np.radians(l)
D=lambda t,n:np.column_stack([np.ones_like(t)]+[f(k*t) for k in range(1,n+1) for f in(np.cos,np.sin)])
def fit(y,n):
 A=D(th,n)*np.sqrt(w)[:,None];return np.linalg.lstsq(A,y*np.sqrt(w),rcond=None)[0]
p1=fit(r,1);r2=r-(p1[1]*np.cos(th)+p1[2]*np.sin(th))
g=np.linspace(0,360,721);gt=np.radians(g);cols=['tab:blue','tab:red','tab:green','tab:orange']
fig,axs=plt.subplots(2,1,figsize=(13,14),sharex=True)
for ax,y,ttl in((axs[0],r,'A: residuals (axisym model removed), n=1–4 components'),(axs[1],r2,'B: after removing n=1, rerun: n=1–4 components')):
 p=fit(y,6);ax.errorbar(l,y,yerr=s,fmt='o',color='k',ms=4,alpha=.6,capsize=2,label='residual')
 for n in range(1,5):
  a,b=p[2*n-1],p[2*n];C=np.hypot(a,b);ph=np.degrees(np.arctan2(b,a))%360/n
  ax.plot(g,a*np.cos(n*gt)+b*np.sin(n*gt),color=cols[n-1],lw=2.5,label=f'n={n}: C={C:.3f}, ℓ_n={ph:.0f}°')
 ax.plot(g,D(gt,4)@p[:9],'k',lw=2,ls='--',label='sum n≤4 + c₀')
 ax.axhline(0,color='grey');ax.set_title(ttl);ax.set_ylabel('residual');ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.07),ncol=3,fontsize=13,frameon=False)
axs[1].set_xlabel('Galactic longitude ℓ [deg]');axs[1].set_xlim(0,360)
fig.tight_layout(h_pad=3);fig.savefig('harmonic_components_before_after_n1.png',dpi=130)
