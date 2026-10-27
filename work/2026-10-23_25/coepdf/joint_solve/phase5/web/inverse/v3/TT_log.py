# Log-l TT of the v3 inverse web (FROZEN field, see tt_inverse3.py) vs Planck 2018 TT with asymmetric error bars, plus v2 seed-1 end and the
# earlier shell-concentrated web. Usage: python3 TT_log.py [dir/]. Fails on NaN.
import sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
W='/workspace/coepdf/joint_solve/phase5/web/inverse/v3/'+(sys.argv[1] if len(sys.argv)>1 else '');I='/workspace/coepdf/joint_solve/phase5/web/inverse/'
z=np.load(W+'tt_inverse3.npz');Pf=z['Pf'];r=json.load(open(W+'tt_inverse3.json'))['results'];v2=np.load(I+'tt_inverse.npz')['s1_end'];v2r=json.load(open(I+'tt_inverse.json'))['results']['s1_end']
b=np.load(I+'../ablation/ablation_spectra.npz')['base'];l=np.arange(len(b))
fig,ax=plt.subplots(figsize=(11,6));ax.errorbar(Pf[:,0],Pf[:,1],yerr=[Pf[:,2],Pf[:,3]],fmt='.',color='0.6',ms=2,elinewidth=0.5,label='Planck 2018 TT',zorder=1)
ax.plot(l[2:],b[2:],color='r',lw=0.8,alpha=0.6,label='earlier shell-concentrated web (bump at l~24)')
ax.plot(l[2:],v2[2:],color='m',lw=0.8,alpha=0.6,label=f"v2 seed 1 (positions only, l<=60): chi2={v2r['chi2_total']:.3g}")
ks=[k for k in z.files if k.endswith('_end')]
for k in ks:
    if not np.all(np.isfinite(z[k])):raise SystemExit('NONFINITE '+k)
    ax.plot(l[2:],z[k][2:],lw=1.2,label=f"v3 {k[:-4]}: chi2={r[k]['chi2_total']:.4g} (l2-30 {r[k]['chi2_2_30']:.4g}, l30-2500 {r[k]['chi2_30_2500']:.4g})")
top=max([3e4]+[2*float(z[k][2:129].max()) for k in ks]);ax.set_xscale('log');ax.set_yscale('log');ax.set_ylim(10,top);ax.set_xlim(2,2500)
ax.axvline(130,color='gray',ls=':',lw=0.8);ax.set_xlabel('multipole l (dotted: l=130 fit limit)');ax.set_ylabel('D_l [uK^2]')
ax.set_title('TT: v3 inverse web (position + mass free, fit to l<=130), FROZEN field, not refitted');ax.legend(fontsize=7);fig.tight_layout();fig.savefig(W+'TT_inverse3_log.png',dpi=130)
