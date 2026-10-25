import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
J=json.load(open('shrink_results.json'));R=J['runs'];nl=J['n_law']
fig,ax=plt.subplots(figsize=(16,10))
for key,lab,c in (('b_old','with Migkas-fitted mass, scaled like old setup (mean 1+z ≈ 2)','C1'),('b_1100','with Migkas-fitted mass, mean 1+z ≈ 1100','C3'),('none_old','no mass, mean 1+z ≈ 2','C0'),('none_1100','no mass, mean 1+z ≈ 1100 (steepest run did not resolve)','C2')):
    ks=[k for k in ('0','0.3','1','3') if f'{key}_k{k}' in R];n=[R[f'{key}_k{k}']['slope_n'] for k in ks];y=[R[f'{key}_k{k}']['xPlanck'] for k in ks]
    ax.plot(n,y,'o-',color=c,lw=3,ms=10,label=f"{lab}: CatWISE-locked {R[key+'_k0']['xPlanck']:.2f}× Planck")
    for ni,yi,k in zip(n,y,ks): ax.annotate(f"1+z {R[f'{key}_k{k}']['mean_1pz']:.0f}",(ni,yi),xytext=(6,-16 if c in('C0','C2') else 6),textcoords='offset points',fontsize=11,color=c)
ax.axhline(1,color='k',ls='--',lw=2);ax.text(2.2,4,'Planck level (1×)',fontsize=14)
ax.axvline(nl,color='m',lw=2);ax.text(nl+.03,200,f'CatWISE 3% lock: slope {nl:.3f}',color='m',fontsize=14)
ax.text(2.55,120,'slope 3 = density falls as 1/r³:\nradius stops shrinking\n(old fixed-radius case)',fontsize=13)
ax.set_xlabel('outward density slope today (−d ln ρ / d ln r at our distance): sets how fast the recombination radius sweeps inward',fontsize=15)
ax.set_ylabel('CMB dipole (× Planck 3362 μK)',fontsize=16);ax.tick_params(labelsize=14);ax.grid(alpha=.3)
ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),fontsize=13,frameon=False)
fig.savefig('dipole_vs_shrink.png',dpi=90,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(16,10))
for f,lab,c in (('b_old_k0','with mass, CatWISE slope','C1'),('none_old_k0','no mass, CatWISE slope','C0'),('none_old_k1','no mass, slope 1.04 (for comparison)','C2')):
    Z=np.load(f+'.npz');t0=float(Z['t0']);i=np.argsort(Z['te']);ax.plot(Z['te'][i]/t0,Z['Rr'][i],color=c,lw=3,label=f'{lab}: radius sweeps inward at {-R[f]["dRrec_dt"]:.0f}× light speed')
ax.axhline(.5,color='k',ls='--',lw=2,label='old fixed recombination radius 0.5 D')
ax.set_xlabel('time / age today (each curve its own age)',fontsize=16);ax.set_ylabel('recombination radius (units of our distance D)',fontsize=16);ax.tick_params(labelsize=14);ax.grid(alpha=.3)
ax.set_ylim(0,1.0);ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),fontsize=13,frameon=False)
fig.savefig('Rrec_vs_time.png',dpi=90,bbox_inches='tight')
