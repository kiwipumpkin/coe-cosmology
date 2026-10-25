import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
R=json.load(open('floor_results.json'));rs=[0.5,0.2,0.05,0.01,0.002,0.000909]
fig,ax=plt.subplots(figsize=(16,10))
for m,lab,c in (('b','Migkas-fitted mass (b), simultaneous emission','C1'),('a','smaller equal mass (a), simultaneous emission','C0')):
    y=[R[f'{m}_{r}']['xPlanck'] for r in rs];ax.plot(rs,y,'o-',color=c,lw=3,ms=10,label=f'{lab}: {y[0]:.2f}× at 0.5 D → {y[-1]:.2f}× at D/1100')
    for r,v in zip(rs,y):ax.annotate(f"1+z {R[f'{m}_{r}']['mean_1pz']:.0f}",(r,v),xytext=(5,8),textcoords='offset points',fontsize=12,color=c)
ax.plot([0.5,0.000909],[R['n_0.5']['xPlanck'],R['n_0.000909']['xPlanck']],'s--',color='C2',lw=3,ms=10,label='no mass, simultaneous emission: 0.00× (exactly flat sky)')
ax.plot([0.000909],[R['shrink_b_1100']['xPlanck']],'*',color='C3',ms=26,label=f"shrinking radius, mass b, CatWISE slope (earlier run): {R['shrink_b_1100']['xPlanck']:.2f}× — same as the floor")
ax.axhline(1,color='k',ls='--',lw=2);ax.text(0.0012,1.1,'Planck dipole = our 370 km/s motion (1×)',fontsize=14)
ax.axhline(.1,color='gray',ls=':',lw=2);ax.text(0.0012,.15,'0.1× target',fontsize=13,color='gray')
ax.set_xscale('log');ax.set_xlabel('recombination radius of our own shell / D (smaller = earlier emission, higher 1+z)',fontsize=16);ax.set_ylabel('model-only CMB dipole (× Planck 3362 μK)',fontsize=16)
ax.tick_params(labelsize=14);ax.grid(alpha=.3,which='both');ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),fontsize=13,frameon=False)
fig.savefig('floor_vs_Rrec.png',dpi=90,bbox_inches='tight')
