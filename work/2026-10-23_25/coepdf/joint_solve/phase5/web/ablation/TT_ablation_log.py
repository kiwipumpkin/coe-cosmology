# Overlay (log l) of web TT ablation variants vs Planck 2018 TT; data from ablation_spectra.npz (made by ablation.py).
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase5/web/ablation/';z=np.load(O+'ablation_spectra.npz');Pf=z['Pf'];l=np.arange(len(z['base']))
for k in ['base','noweb','noISW','nonear','obsonly']:
    if not np.all(np.isfinite(z[k])):raise SystemExit('NONFINITE '+k)
fig,ax=plt.subplots(figsize=(11,6));ax.errorbar(Pf[:,0],Pf[:,1],yerr=[Pf[:,2],Pf[:,3]],fmt='.',color='0.6',ms=2,elinewidth=0.4,label='Planck 2018 TT',zorder=1)
st=dict(base=('k','-','baseline web (RA eps=0.01, M1e22)'),noweb=('r','-','(1) no web'),noISW=('b','--','(2) web, ISW off'),nonear=('g',':','(3) web, nearest 10% bodies removed'),obsonly=('m','-.','diag: Dfix = target observed sky'))
for k,(c,ls,lab) in st.items():ax.plot(l[2:],z[k][2:],color=c,ls=ls,lw=1.2,label=lab,zorder=3)
ax.set_xscale('log');ax.set_yscale('symlog',linthresh=100);ax.set_xlim(2,2500);ax.set_ylim(0,1e6)
ax.axvspan(15,35,color='y',alpha=0.15);ax.set_xlabel(r'$\ell$');ax.set_ylabel(r'$D_\ell$ [$\mu$K$^2$]');ax.legend(fontsize=8,loc='upper right')
ax.set_title('Web TT ablations (frozen field); shaded l=15-35 bump');fig.tight_layout();fig.savefig(O+'TT_ablation_log.png',dpi=130)
