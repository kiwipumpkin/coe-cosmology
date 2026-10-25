import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
S=json.load(open('m345_summary.json'));fig,ax=plt.subplots(figsize=(16,10))
for f,col,lab in (('mass','#cb181d','masses only'),('void','#3182bd','voids only'),('both','#6a3d9a','masses + voids')):
    c=np.array(S[f]['curve']);a=S[f]['after_1e5']
    ax.plot(c[:,0],c[:,1],'o-',color=col,lw=3,ms=10,label=f"{lab}: 1× at N = {S[f]['first_N_1x']}, best {S[f]['best']:.3f}×; side effect: l = 6–12 rises 0.55× → {a['l6_12']:.0f}×")
ax.axhline(1,color='k',ls='--',lw=2.5,label='1× target (10⁻⁵ = 27 μK RMS)')
ax.set_yscale('log');ax.set_xlabel('number of bodies N',fontsize=17);ax.set_ylabel('RMS of l = 3, 4, 5 (× 10⁻⁵)',fontsize=17);ax.tick_params(labelsize=14);ax.grid(alpha=.3,which='both')
ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),fontsize=13,frameon=False,title=f"dipole removed by motion ({S['start']['v_kms']:.0f} km/s); start l=3,4,5 = {S['start_1e5']['l3']:.1f}, {S['start_1e5']['l4']:.2f}, {S['start_1e5']['l5']:.2f} × 10⁻⁵; quadrupole 366× left alone",title_fontsize=14)
fig.text(.5,-.04,'Bodies cut the recombination surface at 14.3 Gpc; dT/T = ln(1+δ)/3 over the patch. Limits from Chris\'s seed-3 laws. Greedy addition + joint refit; dipole re-absorbed by motion after.',ha='center',fontsize=12)
fig.savefig('leftover_m345_vs_N.png',dpi=85,bbox_inches='tight')
