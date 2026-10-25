import json,glob,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':22})
S=[json.load(open(f)) for f in glob.glob('gr2/v4/shadow/shM_*.json')]+[json.load(open('gr2/v4/shadow/sh_eta_fwd_gtc_1.json'))]
S=[d for d in S if 'M' in d]; S.sort(key=lambda d:d['M'])
M=np.array([d['M'] for d in S]); W=np.array([d['shadow_width_deg'] for d in S]); C=np.array([d['chi2'] for d in S]); U=np.array([d['Uo_c'] for d in S])
fig,ax=plt.subplots(2,1,figsize=(19.2,16),gridspec_kw=dict(height_ratios=[3,1.3]))
ax[0].plot(M,W,'o-',c='C0',lw=4,ms=14,label='eta fwd U>c: M varied, other 5 params refit to Migkas (K held 8.46 K)')
for k,(m,w,c,u) in enumerate(zip(M,W,C,U)): ax[0].annotate(f'χ²={c:.1f}, U/c={u:.2g}',(m,w),textcoords='offset points',xytext=(-330,-8-0*k),fontsize=14,arrowprops=dict(arrowstyle='-',lw=0.8))
mk=dict(eta_fwd_gtc_1='*',pl_back_gtc_3='s',pl_fwd_le_1='D',soft_fwd_gtc_0='^',exp_back_le_0='v')
for t,s in mk.items():
    d=json.load(open(f'gr2/v4/shadow/sh_{t}.json')); ax[0].plot(d['M'],max(d['shadow_width_deg'],0.03),s,ms=22,mec='k',label=f"{t}: {d['shadow_width_deg']:.2f}°")
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_ylabel('CoE shadow full width (deg)'); ax[0].set_xlabel('fitted mass M (code units)'); ax[0].grid(which='both'); ax[0].legend(fontsize=15,loc='lower right')
ax[0].set_title('CoE shadow: directions whose past ray ends in the mass without reaching 3000 K (2° grid, edges bisected to 0.06°)',fontsize=20)
d=json.load(open('gr2/v4/shadow/sh_eta_fwd_gtc_1.json')); g=[x for x in d['grid'] if x['rec']]; e=np.array([x['ell'] for x in g]); p=np.array([x['phi_e_deg'] for x in g]); T=np.array([x['T_K'] for x in g]); j=np.where(np.diff(e)>3)[0]+1; e=np.insert(e.astype(float),j,np.nan); p=np.insert(p,j,np.nan); T=np.insert(T,j,np.nan)
a2=ax[1]; a2.plot(e,np.abs(p),'o-',lw=3,label='|emitter azimuth| (deg): >180 = image from behind the mass'); a2.axvspan(*d['edges'],color='0.5',alpha=0.4); a2.axhline(180,c='k',ls='--')
a3=a2.twinx(); a3.plot(e,T,'s-',c='C3',lw=3,label='T_obs (K)'); a3.set_ylabel('T_obs (K)',color='C3')
a2.set_xlabel('model angle θ (deg; CoE at 180)  — eta_fwd_gtc_1'); a2.set_ylabel('|φ_e| (deg)'); a2.legend(fontsize=15,loc='upper left'); a2.grid()
fig.tight_layout(); fig.savefig('gr2/v4/shadow/shadow_vs_M.png',dpi=100)
ok=(C<100)&(U>1.01); sl=np.polyfit(np.log10(M[ok]),np.log10(W[ok]),1)[0] if ok.sum()>1 else np.nan
print(list(zip(M,W,C,U)),'slope',sl)
