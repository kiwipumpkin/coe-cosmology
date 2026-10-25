import sys,json; sys.path.insert(0,'/workspace/cosmology-course/notes/paper_repro/grfit'); from core import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
p=json.load(open('/workspace/cosmology-course/notes/paper_repro/grfit/fit_newt_plain.json'))['rows'][0]
af=nsky(p['D'],p['mcc'],p['kin'],p['kout']); ML=TH1[IMIG]
fig,ax=plt.subplots(figsize=(9,5.5))
ax.errorbar(ML,Y,yerr=[MLO,MHI],fmt='o',ms=3,color='k',lw=.8,label='Migkas et al. 2020 Fig. 6, A/A_all (71 bins)')
for l0,c,ls,nm in [(p['l0'],'C0','-','paper Newtonian best fit, l0=%.3f°'%p['l0']),(15.166,'C3','--','assumed anti-CoE l=15.166° (b=+30° does not enter)'),(195.166,'C2',':','flipped: CoE at l=15.166°')]:
    m=shifted(af,l0); print(nm,chi(m),r2(m)); ax.plot(TH1,m,color=c,ls=ls,lw=1.6,label=f"{nm}: χ²={chi(m):.2f}, R²={r2(m):.3f}")
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel('Galactic longitude l (deg); model angle from anti-CoE = l − l0'); ax.set_ylabel('A / A_all')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=8,frameon=False); fig.tight_layout(); fig.savefig('../migkas_A_vs_angle_coe.png',dpi=150)
print(p['D'],p['chi2'],p.get('R2'))
