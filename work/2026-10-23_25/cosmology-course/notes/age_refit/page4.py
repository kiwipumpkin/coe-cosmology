import json,sys,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':20})
d=json.load(open(sys.argv[1])); out=sys.argv[2]
fig,ax=plt.subplots(3,1,figsize=(19.2,21.6))
md=d['migkas_data']; ax[0].errorbar(md['l'],md['y'],yerr=md['err'],fmt='o',c='k',ms=7,label='Migkas A/A_all (71 bins)')
ax[0].plot(np.arange(360),d['migkas_model'],c='C3',lw=4,label=f"exact-GR model  χ²={d['chi2']:.2f}  σ={d['sigma']:.1f}  l0={d['l0']:.0f}°")
ax[0].set_ylabel('A / A_all'); ax[0].legend(fontsize=20)
sp=d.get('score_partial'); 
if sp:
    gl=np.array(sp['gal_l']); i=np.argsort(gl)
    ax[1].plot(gl[i],100*np.array(sp['dTT'])[i],'o-',c='C0',lw=4,ms=9,label=f"model CMB dT/T (T0 solved to 2.725 K mean)")
    import importlib.util as u; s=u.spec_from_file_location('cc','/workspace/clump_scan/cmb-cancel.py'); m=u.module_from_spec(s); s.loader.exec_module(m)
    ax[1].plot(np.arange(360),-100*(m.G_FULL-1),c='k',lw=2.5,ls='--',label='needed for flat CMB: −(G−1)')
    ax[1].set_ylabel('dT/T  (%)'); ax[1].legend(fontsize=18)
    ax[2].plot(gl[i],100*np.array(sp['resid'])[i],'s-',c='C2',lw=4,ms=9,label=f"residual G·S, m=1 removed: rms={sp['res_1e5']/1e3:.1f}×10⁻² ; cancellation {sp['cancel']:.0f}%")
    ax[2].axhline(0,c='k'); ax[2].set_ylabel('residual (%)'); ax[2].legend(fontsize=18)
for r in d.get('rows',[]):
    if r['status']!='rec':
        for a in ax[1:]: a.axvspan((r['ell']+d['l0']-3)%360,(r['ell']+d['l0']+3)%360,color='0.5',alpha=0.4)
for a in ax: a.set_xlim(0,360); a.grid(); a.tick_params(labelsize=18)
ax[2].set_xlabel('galactic longitude ℓ (deg)   [CoE at l0+180; grey = no recombined emitter on the ray (core shadow)]')
T=d.get('T_obs_range_K',[np.nan,np.nan])
fig.suptitle(f"{d['outtag']}:  D={d['D_Gpc']:.3g} Gpc   U_o/c={d['Uo_c']:.3g}   M={d['M']:.3g}\nK(local radiation T today)={d.get('K_local_today_K',np.nan):.4g} K   <T_obs>={d.get('T_obs_mean_K',np.nan):.4g} K   range {T[0]:.3g}–{T[1]:.3g} K   coverage {d.get('coverage')}/{d['n_dir']}",fontsize=22)
fig.tight_layout(); fig.savefig(out,dpi=100)
