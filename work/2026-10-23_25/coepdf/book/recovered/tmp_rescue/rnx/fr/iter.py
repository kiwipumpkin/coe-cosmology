import numpy as np,json
from scipy.stats import pearsonr,chi2 as C2
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('residuals.csv',delimiter=',',names=True)
l,r,s=d['ell_deg'],d['residual'],d['sigma'];w=1/s**2;th=np.radians(l);N=len(l)
D=lambda t,n:np.column_stack([np.ones_like(t)]+[f(k*t) for k in range(1,n+1) for f in(np.cos,np.sin)])
def fit(y,X):
 A=X*np.sqrt(w)[:,None];cov=np.linalg.inv(A.T@A);p=cov@(A.T@(y*np.sqrt(w)));return p,cov,float(np.sum(w*(y-X@p)**2))
p1,_,_=fit(r,D(th,1)); n1=p1[1]*np.cos(th)+p1[2]*np.sin(th)
r2=r-n1; out={'n1_removed':dict(a1=p1[1],b1=p1[2],C1=np.hypot(p1[1],p1[2]),peak=np.degrees(np.arctan2(p1[2],p1[1]))%360),'chi2_after_n1':float(np.sum(w*r2**2))}
prev=fit(r2,D(th,0))[2];out['chi2_offset']=prev;tab=[]
pF,cv,_=fit(r2,D(th,6))
for n in range(1,7):
 p,_,c=fit(r2,D(th,n));a,b=pF[2*n-1],pF[2*n];C=np.hypot(a,b)
 sC=np.sqrt(a*a*cv[2*n-1,2*n-1]+b*b*cv[2*n,2*n]+2*a*b*cv[2*n-1,2*n])/C
 tab.append(dict(n=n,C=C,sC=sC,peak=(np.degrees(np.arctan2(b,a))%360)/n,chi2=c,dchi2=prev-c,p=float(C2.sf(prev-c,2)),r2=pearsonr(r2,D(th,n)@p)[0]**2));prev=c
out['table']=tab;json.dump(out,open('iter_after_n1.json','w'),indent=1,default=float)
print(json.dumps(out,indent=1,default=float))
g=np.linspace(0,360,721);gt=np.radians(g)
plt.rcParams.update({'font.size':16})
fig,ax=plt.subplots(figsize=(13,7.5))
ax.errorbar(l,r2,yerr=s,fmt='o',color='k',ms=5,capsize=2,label='residual after removing n=1')
for n,c in((2,'tab:red'),(3,'tab:blue')):
 p,_,cc=fit(r2,D(th,n));ax.plot(g,D(gt,n)@p,color=c,lw=2.5,label=f'new Fourier fit n≤{n} (χ²={cc:.1f})')
ax.axhline(0,color='grey');ax.set_xlim(0,360);ax.set_xlabel('Galactic longitude ℓ [deg]');ax.set_ylabel('r(ℓ) − C₁cos(ℓ−ℓ₁)')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.13),ncol=3,fontsize=13,frameon=False);fig.tight_layout();fig.savefig('residual_after_n1_fourier.png',dpi=150)
