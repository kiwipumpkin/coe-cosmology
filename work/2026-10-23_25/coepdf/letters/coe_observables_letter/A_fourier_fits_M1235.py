# Naive truncated Fourier fit of Migkas et al. 2020 Fig.6 A/A_all vs Galactic longitude (71 points, 5 deg step).
# Weights 1/sigma^2 with sigma = mean of asymmetric errors (el+eh)/2. Model: a0 + sum_{m=1..M} [a_m cos(m l) + b_m sin(m l)].
# NOTE: the 71 points come from overlapping 90-deg-wide regions (Migkas Sect.5.2), so they are not independent;
# chi2/dof and AIC/BIC below treat them as independent, as a naive fit does.
import numpy as np,matplotlib,json;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high']
s=(el+eh)/2;w=1/s**2;L=np.radians(l);n=len(l)
def X(M,LL):return np.stack([np.ones_like(LL)]+sum([[np.cos(k*LL),np.sin(k*LL)] for k in range(1,M+1)],[]),1)
def fit(M):
  A=X(M,L);co=np.linalg.solve(A.T@(A*w[:,None]),A.T@(w*y));r=y-A@co;chi=float((w*r*r).sum())
  ybar=(w*y).sum()/w.sum();R2=1-chi/float((w*(y-ybar)**2).sum());k=2*M+1
  return co,chi,R2,k,chi+2*k,chi+k*np.log(n)
rows=[];res={}
for M in range(1,9):
  co,chi,R2,k,aic,bic=fit(M);res[M]=co;rows.append(dict(M=M,k=k,chi2=chi,dof=n-k,chi2_dof=chi/(n-k),R2w=R2,AIC=aic,BIC=bic))
a0=min(r['AIC'] for r in rows);b0=min(r['BIC'] for r in rows)
for r in rows:r['dAIC']=r['AIC']-a0;r['dBIC']=r['BIC']-b0;print({k:(round(v,3) if isinstance(v,float) else v) for k,v in r.items()})
co=res[2];coef={'a0':float(co[0])}
for m in (1,2):
  a,b=co[2*m-1],co[2*m];amp=np.hypot(a,b);ph=np.degrees(np.arctan2(b,a))%360
  coef['m%d'%m]={'amp':float(amp),'phase_deg_of_max':float(ph/m),'a':float(a),'b':float(b)}
print('M=2',coef);json.dump({'table':rows,'M2_coeffs':coef},open('A_fourier_naive_fit.json','w'),indent=1)
g=np.linspace(0,360,721);G=np.radians(g);fig,ax=plt.subplots(2,1,figsize=(11,10))
for a,Ms in zip(ax,[[2],[1,2,3,5]]):
  a.errorbar(l,y,yerr=[el,eh],fmt='o',color='gray',ms=3,capsize=2,label='Migkas 2020 Fig.6 (digitised)')
  for M,c in zip(Ms,['crimson','royalblue','green','purple'] if len(Ms)>1 else ['crimson']):
    r=rows[M-1];a.plot(g,X(M,G)@res[M],color=c,lw=1.6,label='M=%d: R$^2$=%.3f, $\\chi^2$=%.2f, $\\chi^2$/dof=%.3f'%(M,r['R2w'],r['chi2'],r['chi2_dof']))
  a.set(xlabel='Galactic longitude l [deg]',ylabel='A/A$_{all}$',xlim=(0,360));a.grid(alpha=.3);a.legend(fontsize=8)
c1,c2=coef['m1'],coef['m2']
ax[0].set_title('Weighted Fourier fit m=0,1,2: a0=%.4f; m=1 amp %.4f max at l=%.1f$^\\circ$; m=2 amp %.4f max at l=%.1f$^\\circ$ (+180)'%(coef['a0'],c1['amp'],c1['phase_deg_of_max'],c2['amp'],c2['phase_deg_of_max']),fontsize=9)
ax[1].set_title('Fits with M = 1, 2, 3, 5')
fig.text(.5,.004,'Weights 1/$\\sigma^2$, $\\sigma$=(el+eh)/2. Points come from overlapping 90$^\\circ$ regions, so not independent; $\\chi^2$/dof treats them as independent.',ha='center',fontsize=8)
plt.tight_layout(rect=[0,.015,1,1]);plt.savefig('A_fourier_naive_fit.png',dpi=140)
