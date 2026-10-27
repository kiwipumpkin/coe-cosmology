# Log-l TT of the inverse-problem web (FROZEN field, see tt_inverse.py) vs Planck 2018 TT with asymmetric error bars; also old web baseline.
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
W='/workspace/coepdf/joint_solve/phase5/web/inverse/';z=np.load(W+'tt_inverse.npz');Pf=z['Pf'];r=json.load(open(W+'tt_inverse.json'))['results']
b=np.load(W+'../ablation/ablation_spectra.npz')['base'];l=np.arange(len(b))
fig,ax=plt.subplots(figsize=(11,6));ax.errorbar(Pf[:,0],Pf[:,1],yerr=[Pf[:,2],Pf[:,3]],fmt='.',color='0.6',ms=2,elinewidth=0.5,label='Planck 2018 TT',zorder=1)
ax.plot(l[2:],b[2:],color='r',lw=1,alpha=0.7,label='old web compensation (bump at l~24)')
for k in [k for k in z.files if k.endswith('_end')]:
    if not np.all(np.isfinite(z[k])):raise SystemExit('NONFINITE '+k)
    ax.plot(l[2:],z[k][2:],lw=1,label=f"inverse web {k[:-4]}: chi2={r[k]['chi2_total']:.0f} (l2-30 {r[k]['chi2_2_30']:.0f})")
ax.set_xscale('log');ax.set_yscale('log');ax.set_ylim(10,max(3e4,2*max(float(z[k][2:129].max()) for k in z.files if k.endswith('_end'))));ax.set_xlim(2,2500);ax.set_xlabel('multipole l');ax.set_ylabel('D_l [uK^2]')
ax.set_title('TT: uniform random web + minimal displacements (FROZEN field, not refitted)');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(W+'TT_inverse_log.png',dpi=130)
