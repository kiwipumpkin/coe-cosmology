# PROVENANCE COPY of /workspace/coepdf/d5_mod_earth.py (original generating script, rev6 figure; run in its original directory).
# Output figure copied here as L16_d5_mod.png. Not re-run for rev7: see Appendix B.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11})
exec(open('/workspace/coepdf/scripts/shp.py').read())
r0,dv=42.0,14.3005
m=lambda r:shp(r/42)
s=(m(r0+1e-3)-m(r0-1e-3))/2e-3
fig,ax=plt.subplots(figsize=(13,4.6));r=np.linspace(0,80,600)
ax.axvspan(r0-dv,r0+dv,color='gold',alpha=0.3,label=f'observable range $r_\\oplus\\pm${dv:.1f} Gpc ({r0-dv:.1f}–{r0+dv:.1f} Gpc)')
re=r0+dv;mf=lambda r:np.where(r>re,m(re),m(r))
ax.axvspan(re,80,color='0.85',alpha=0.5,label='beyond observable edge: unconstrained')
se=(m(re)-m(re-1e-3))/1e-3;x=np.linspace(0,80-re,200);ys=[]
for i,k in enumerate([2e-5,6e-5,1.5e-4,3e-4]):
  y=m(re)+se*x-k*x**2;ys+=list(y);ax.plot(re+x,y,ls='--',lw=1.5,color=plt.cm.Blues(0.45+0.15*i),label='illustrative continuations (unconstrained)' if i==0 else None)
ri=r[r<=re];ax.plot(ri,m(ri),'b',lw=2,label=r'$m(r)=\exp\sum_k c_k(\min(r/42,R_b)-1)^{k+1}$, $c$=[−0.041, −0.035, 0.1353, 0], $R_b$=2.011')
rt=np.linspace(r0-1.6*dv,r0+1.6*dv,10);ax.plot(rt,m(r0)+s*(rt-r0),'r--',lw=1.6,label=f'tangent at Earth: slope {s:.2e} per Gpc')
ax.plot(r0,m(r0),'ko',ms=9,label=f'Earth, r = {r0:.0f} Gpc, m = {m(r0):.3f}')
ax.set_xlim(0,80);ax.set_ylim(0.97*min(min(m(ri)),min(ys)),1.03*max(max(m(ri)),max(ys)));ax.set_xlabel('distance from CoE axis point [Gpc]');ax.set_ylabel('m (before normalising)');ax.set_facecolor('w')
fig.legend(loc='lower center',ncol=2,frameon=False,fontsize=9);fig.tight_layout(rect=(0,0.16,1,1));fig.savefig('fig/d5_mod.png',dpi=150,facecolor='w')
print(s,m(r0),m(r0-dv),m(r0+dv))
