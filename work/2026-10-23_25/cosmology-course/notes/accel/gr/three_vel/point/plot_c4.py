exec(open('point_c4.py').read().split("# multistart")[0])
import json
from scipy.stats import chi2 as C2,norm
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
P=json.load(open('prof_c4.json'));xc=np.array(P['5'])if False else np.array(P['8']['x'])
xu=np.array([json.load(open('point_fit.json'))['best'][n] for n in ['core_frac','log10_uD','l0','beta_f','v_pec_kms','l_pec','sin_b_pec']])
wm=np.average(Y,weights=1/S**2);cf=float(np.sum(((Y-wm)/S)**2))
o=np.argsort(L);fig,ax=plt.subplots(figsize=(13,8))
ax.errorbar(L,Y,yerr=S,fmt='o',color='k',ms=4,alpha=.6,label=f'Migkas A/A_all data ({len(Y)} bins)')
for x,col,lab in ((xu,'C0','unconstrained best fit'),(xc,'C3','constrained (speed ≥ 4 c), best of profile; R² flat from 4 to 30 c')):
    m=model(x);c=float(np.sum(((m-Y)/S)**2));R2=1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2);sg=norm.isf(C2.sf(cf-c,6)/2)
    print(lab,10**x[1],x[0],c,R2,sg)
    ax.plot(L[o],m[o],'-',color=col,lw=2.5,label=f"{lab}: speed {10**x[1]:.3f} c, core mass {x[0]:.4f} of critical, R² {R2:.3f}, χ² {c:.2f}, {sg:.2f}σ vs flat")
ax.grid(alpha=.3);ax.set_xlabel('Galactic longitude of cluster bin (deg)');ax.set_ylabel('A / A_all')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
fig.savefig('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/migkas_c4_fit.png',dpi=110,bbox_inches='tight')
