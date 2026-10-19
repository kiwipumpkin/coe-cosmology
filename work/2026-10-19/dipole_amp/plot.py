import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from an import dip
S=[0,.25,.5,.75,.9,.95,.99,1];R={}
for s in S: R[s]=dip('chk.npz' if s==0 else f's_{s}.npz')
x=np.array(S);y=np.array([R[s]['xPlanck'] for s in S]);spr=np.array([(R[s]['te_max']-R[s]['te_min'])/(0.5*(R[s]['te_max']+R[s]['te_min']))*100 for s in S])
k,b=np.polyfit(1-x,y,1);s1=1-(1-b)/k
out=dict(simultaneity={str(s):dict(R[s],emission_time_spread_pct=float(p)) for s,p in zip(S,spr)},linear_fit=dict(slope_per_unit_1_minus_s=k,intercept=b,s_for_1x_Planck=s1),
 speed={u:dip(f'u_{u}.npz') for u in [1,2,8,16,30]},rrec={r:dip(f'r_{r}.npz') for r in [0.3,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001]},
 reference={k2:dip(f) for k2,f in [('no mass own','../z4_nomass.npz'),('equal mass a own','../distributed/z4_dist_a_own.npz'),('b own','../distributed/z4_dist_b_own.npz'),('b sim','../distributed/z4_dist_b_sim.npz'),('a sim','../distributed/z4_dist_a_sim.npz')]})
json.dump(out,open('dipole_amp.json','w'),indent=1,default=float)
fig,ax=plt.subplots(figsize=(16,10))
ax.plot(x,y,'o-',color='C1',lw=3,ms=10,label=f'computed (own crossing, distributed mass b, 4c, recombination radius 0.5 D): dipole ≈ {b:.2f} + {k:.1f} × (1 − simultaneity)')
ax.axhline(1,color='k',ls='--',lw=2,label='Planck dipole level (3362 μK)')
ax.axvline(s1,color='C2',lw=2,label=f'1× Planck reached at simultaneity {s1:.4f} (emission-time spread {np.interp(s1,x,spr):.2f}% vs {spr[0]:.1f}% fully sequential)')
for s in (0,.5,.9,1): ax.annotate(f"{R[s]['xPlanck']:.1f}× Planck\nspread {np.interp(s,x,spr):.1f}%\nMigkas χ² 7.673 (unchanged)",(s,R[s]['xPlanck']),textcoords='offset points',xytext=(10,15),fontsize=13)
ax.set_xlabel('simultaneity: 0 = each emitter at its own crossing of 0.5 D, 1 = all emit at one moment',fontsize=16);ax.set_ylabel('dipole amplitude (× Planck 3362 μK)',fontsize=16)
ax.tick_params(labelsize=14);ax.grid(alpha=.3);ax.set_xlim(-.02,1.05);ax.set_ylim(-3,100)
ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),fontsize=13,frameon=False)
ia=ax.inset_axes([.5,.42,.35,.3]);ia.plot(x,y,'o-',color='C1',lw=2);ia.axhline(1,color='k',ls='--');ia.axvline(s1,color='C2');ia.set_xlim(.94,1.002);ia.set_ylim(0,6);ia.set_title('zoom near full simultaneity',fontsize=12);ia.grid(alpha=.3)
fig.savefig('dipole_vs_simultaneity.png',dpi=90,bbox_inches='tight');print(k,b,s1,spr)
