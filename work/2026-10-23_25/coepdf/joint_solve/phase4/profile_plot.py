# Phase4 (NEW): profile chi2 vs gradient ratio k (=fitted delta / CatWISE 0.031) and how D, d, L, eta, mc, l0 move along it.
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
P=json.load(open('/workspace/coepdf/joint_solve/phase4/profile_k.json'));out={}
fig,ax=plt.subplots(2,3,figsize=(14,7.5));ax=ax.ravel()
for m,lab,c in (('A','d free [2,40] Gpc','C3'),('B','d fixed 14.30 Gpc','C0')):
    R=[r for r in P if r[0]==m and r[2]<800];k=np.array([10**r[1] for r in R]);ch=np.array([r[3]['chi2'] for r in R]);mn=ch.min()
    lk=np.log10(k);g=np.linspace(lk.min(),lk.max(),2001);cg=np.interp(g,lk,ch);lo1=10**g[np.argmax(cg<mn+1)];lo4=10**g[np.argmax(cg<mn+4)]
    out[m]=dict(chi2_min=float(mn),k_at_min=float(k[ch.argmin()]),k_low_1sigma=float(lo1),k_low_2sigma=float(lo4),dchi2_at_k1=float(np.interp(0,lk,ch)-mn),
        rows=[dict(k=float(kk),**{q:R[i][3][q] for q in ('chi2','r2','D','d','L','eta','mc','l0','age_Gyr')}) for i,kk in enumerate(k)])
    ax[0].semilogx(k,ch,'o-',color=c,label=lab)
    for i,q in enumerate(('D','d','L','eta','mc')):ax[i+1].loglog(k,[r[3][q] for r in R],'o-',color=c,label=lab)
ax[0].axvline(1,color='k',ls=':',label='CatWISE (k=1)');ax[0].set_ylabel('Migkas chi2 (A=(H/<H>)^2), n=71');ax[0].set_ylim(8,46)
for a,t in zip(ax,('profile chi2','D [Gpc]','d [Gpc]','L [Gpc]','eta','mc')):a.set_title(t,fontsize=10);a.set_xlabel('k = fitted gradient / CatWISE');a.grid(alpha=.3);a.legend(fontsize=7)
fig.suptitle('NEW phase4: <H0>=68.64 locked, age 13.787 Gyr ±5% penalty, gradient free',fontsize=11);fig.tight_layout()
fig.savefig(__file__.replace('.py','.png'),dpi=130);json.dump(out,open(__file__.replace('.py','.json'),'w'),indent=1)
print({m:{q:v for q,v in o.items() if q!='rows'} for m,o in out.items()})
