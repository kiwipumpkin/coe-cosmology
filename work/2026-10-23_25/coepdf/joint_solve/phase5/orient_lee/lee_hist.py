import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=__file__.rsplit('/',1)[0];d=json.load(open(H+'/lee_results.json'));r=d['real']
fig,ax=plt.subplots(1,2,figsize=(14,5))
for tag,c in [('isotropic','C0'),('fixedCoE_roll','C1')]:
  h=np.array(d['mc'][tag]['score_hist']);ax[0].bar(np.arange(13)+(-.2 if c=='C0' else .2),h/h.sum(),.4,color=c,label=f"{tag}: p(≥{r['score12']})={d['mc'][tag]['p_score>=real']:.4f}")
  o=np.load(H+f'/off_{tag}.npy');ax[1].hist(o,bins=90,range=(0,90),density=True,histtype='step',lw=2,color=c,label=f"{tag}: p(≤real)={d['mc'][tag]['aoe']['phase4b joint D=27.34/R_ls=21.74']['p_aoe<=real']:.3f}")
ax[0].axvline(r['score12'],color='k',ls='--');ax[0].set_xlabel('sign matches out of 12');ax[0].set_ylabel('fraction');ax[0].legend();ax[0].set_yscale('log')
o0=d['mc']['isotropic']['aoe']['phase4b joint D=27.34/R_ls=21.74']['real_offset_deg']
ax[1].axvline(o0,color='k',ls='--',label=f'paper orientation {o0:.1f}°');ax[1].set_xlabel('AoE offset to nearest critical point (deg), D/R=27.34/21.74');ax[1].legend()
fig.suptitle(f"Orientation look-elsewhere MC, N={d['N']} per ensemble");fig.tight_layout();fig.savefig(H+'/lee_hist.png',dpi=110)
