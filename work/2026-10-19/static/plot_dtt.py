import pickle,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B=pickle.load(open('best.pkl','rb'));res=B['res'];PR=B['PR'];ps=B['psis']
mir=lambda a,y:(np.r_[a,360-a[::-1]],np.r_[y,y[::-1]])
# Figure 1: dT/T 0-360
fig,ax=plt.subplots(figsize=(13,8))
s=res['nomass']['stats'];a,y=mir(ps,s['d'])
ax.plot(a,y,'C0',lw=2.5,label=f"static space, no mass (visible patch only, 150°–210°): RMS {s['rms']:.3f}, range {s['rng']:.3f} ({s['rng']/PR:.0f}× Planck ring); after dipole removal range {s['rrng']:.3f} ({s['rrng']/PR:.0f}× Planck)")
pc=np.arange(0,181.)
for t,col,lab in (('nomass','C2','expanding (comoving) space, no mass, full sky'),('mass','C3','expanding (comoving) space, core mass 0.0154 of critical, full sky')):
    c=res['comoving_'+t];a,y=mir(pc,np.array(c['d']));ax.plot(a,y,col,lw=2,label=f"{lab}: RMS {c['rms']:.3f}, range {c['rng']:.3f} ({c['rng']/PR:.0f}× Planck); after dipole removal range {c['rrng']:.3f} ({c['rrng']/PR:.0f}× Planck)")
ax.axvspan(0,150,color='k',alpha=.07);ax.axvspan(210,360,color='k',alpha=.07)
ax.text(75,0.95,'static space: no light reaches us\nfrom these directions',ha='center',va='top',transform=ax.get_xaxis_transform(),fontsize=11)
ax.text(285,0.95,'static space: no light reaches us\nfrom these directions',ha='center',va='top',transform=ax.get_xaxis_transform(),fontsize=11)
ax.axhline(0,color='k',lw=.6);ax.grid(alpha=.3);ax.set_xlim(0,360);ax.set_xticks(range(0,361,30))
ax.set_xlabel('look direction on the sky from the anti-CoE direction (0°/360° = away from CoE, 180° = toward CoE; mirrored)');ax.set_ylabel('dT/T (relative to mean over the directions that receive light)')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=9.5,frameon=False)
ax.set_title('our speed 1.052 c (Migkas best fit). Static space WITH mass is not drawn: its redshift passes through zero and infinity between 159.2° and 161.9°, so dT/T is undefined',fontsize=9,pad=-14)
fig.savefig('static_dTT_360.png',dpi=110,bbox_inches='tight')
# Figure 2: 1+z vs angle, zoom on visible patch
fig,ax=plt.subplots(figsize=(13,8))
for t,col,lab in (('nomass','C0','static, no mass'),('mass','C1','static, core mass 0.0154 of critical')):
    z=np.array([x['opz'] for x in res[t]['r']]);zz=np.where(np.isfinite(z)&(z>0),z,np.nan)
    a,y=mir(ps,zz);ax.plot(a,y,col,lw=2.5,label=f"{lab}: 1+z = {z[0]:.2f} at patch edge (150°), {z[-1]:.3f} at centre (180°)"+("" if t=='nomass' else "; negative between 159.2° and 161.9°, rises to 1488 at 161.95°"))
for t,col,lab in (('nomass','C2','expanding (comoving), no mass'),('mass','C3','expanding (comoving), core mass 0.0154 of critical')):
    c=res['comoving_'+t];o=np.load('../zmig_'+t+'.npz')['o'][:,1];a,y=mir(pc,o)
    ax.plot(a,y,col,lw=2,ls='--',label=f"{lab}: 1+z = {o[150]:.3f} at 150°, {o[180]:.3f} at 180°")
for lo,hi in ((159.16,161.94),(360-161.94,360-159.16)): ax.axvspan(lo,hi,color='m',alpha=.25)
ax.text(180,14.3,'magenta band: static with mass, 1+z ≤ 0 or infinite (source outruns its own light toward us)',ha='center',fontsize=10)
ax.set_xlim(145,215);ax.set_ylim(0,15);ax.grid(alpha=.3);ax.set_xticks(range(145,216,5))
ax.set_xlabel('look direction from the anti-CoE direction (180° = toward CoE); outside 150°–210° no light reaches us in static space');ax.set_ylabel('1 + z (y-axis cut at 15)')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
fig.savefig('static_redshift_vs_angle.png',dpi=110,bbox_inches='tight')
