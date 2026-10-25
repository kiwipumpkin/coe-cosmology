"""usage: cmb2.py <law> [ckjson]  CMB from dt_step circle emitters, 6-deg, K=2.725 (not tuned). checkpoint per ray."""
import sys,json,pickle,os,numpy as np; from lc3 import *
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
law=sys.argv[1]; src=sys.argv[2] if len(sys.argv)>2 else f'res_{law}.json'; J=json.load(open(src))
x=np.array(J['x']); Ri=10**J['lRi'] if 'lRi' in J else 1e-3; l0=J['l0']; af=np.array(J['af']) if 'af' in J else np.array(json.load(open(f'ck_{law}.json'))['final6']['af'])
L,tn=make2(law,x[:len(NAMES[law])],Ri); ts=tstar(L,tn); o=obs(L,tn)
tag=os.path.basename(src).replace('.json','')
if ts is None: json.dump(dict(law=law,src=src,Ri=Ri,result='our element never reaches 3000 K (t* undefined) -> no emission anywhere'),open(f'cmb3_{os.path.basename(src)[:-5]}.json','w')); sys.exit('no t*')
ck=f'cmb3_{tag}.pkl'; done=pickle.load(open(ck,'rb')) if os.path.exists(ck) else {}
ells=np.arange(0,360,6.)
for e in ells:
    if e in done: continue
    d=emitters_rec(L,tn,ts,[e])[0]; d.pop('pts'); done[e]=d; pickle.dump(done,open(ck,'wb')); print(e,d['status'],flush=True)
rows=[done[e] for e in ells]; ok=[q for q in rows if q['status']=='emitter']
out=dict(law=law,src=src,Ri=Ri,t_star_over_tn=ts/tn,U_c=o['U'],n_emitter=len(ok),n_beyond_CoE=sum(q['status']=='beyond_CoE' for q in rows),
  n_no_emitter=sum(q['status'].startswith('stopped') for q in rows),no_emitter_dirs=[q['ell'] for q in rows if q['status'].startswith('stopped')],
  stops={q['ell']:q['stop'] for q in rows if q['status'].startswith('stopped')},n_no_matter=sum(q.get('no_matter',False) for q in ok))
Tt=np.array([q['T_obs'] for q in ok]); out.update(T_range_K=[Tt.min(),Tt.max()],T_mean_K=Tt.mean(),T_mean_vs_2725=Tt.mean()/2.725,
  z_range=[min(q['z'] for q in ok),max(q['z'] for q in ok)],
  te_over_tn_range=[min(q['t_e'] for q in ok)/tn,max(q['t_e'] for q in ok)/tn])
sc=score([dict(q,ell=q['ell']) for q in ok],l0,ells); out.update({k:sc[k] for k in ('dTT_rms_m1rem','target_rms_m1rem','resid_rms','cancel_pct','corr_with_Gm1')})
te={q['ell']:q['t_e'] for q in ok}; far=max(te,key=lambda e:-abs(e-180)); out['order']=dict(te0=te.get(0.0),ell_coe_side=far,te_coe=te[far],zero_first=bool(te.get(0.0,9e9)<te[far]))
json.dump(out,open(f'cmb3_{tag}.json','w'),indent=1,default=float)
fig,ax=plt.subplots(2,1,figsize=(19,11)); ell=np.arange(360.)
ax[0].errorbar(F.IMIG,F.Y,yerr=0.5*(F.MLO+F.MHI),fmt='o',ms=7,c='k',label='Migkas 71 bins'); ax[0].plot(ell,F.shifted(af,l0),lw=4,c='C3',label=f'light-cone exact GR {law}: χ²={J["chi2"]:.1f}, {J["sigma"]:.1f}σ, l0={l0:.0f}°')
ax[0].set_ylabel('A / A_all',fontsize=22); ax[0].legend(fontsize=18)
gl=np.array(sc['gal_l']); i=np.argsort(gl)
ax[1].plot(gl[i],np.array(sc['dTT'])[i],'o-',lw=3,label=f'model dT/T (m=1 removed), rms {sc["dTT_rms_m1rem"]:.3g}')
ax[1].plot(gl[i],np.array(sc['Gm1'])[i],'k-',lw=3,label='G−1 (m=1 removed)')
ax[1].plot(gl[i],np.array(sc['resid'])[i],':',lw=3,c='C2',label=f'residual {sc["resid_rms"]:.3g}, cancel {sc["cancel_pct"]:.0f}%, corr {sc["corr_with_Gm1"]:+.2f}')
ax[1].text(0.01,0.04,f'{out["n_emitter"]}/60 directions have emitters; {out["n_no_emitter"]} no emitter (horizon/stopped), {out["n_beyond_CoE"]} beyond CoE',transform=ax[1].transAxes,fontsize=16)
ax[1].set_xlabel('galactic longitude ℓ (deg)',fontsize=22); ax[1].set_ylabel('dT/T',fontsize=22); ax[1].legend(fontsize=16)
for a in ax: a.set_xlim(0,360); a.tick_params(labelsize=18); a.grid()
fig.suptitle(f'{law}: U/c={o["U"]:.2f}, T_obs {Tt.min():.2f}–{Tt.max():.2f} K  mean {Tt.mean():.3f} K vs 2.725, T_obs=3000/(1+z), z_e {out["z_range"][0]:.2f}–{out["z_range"][1]:.2f}',fontsize=18)
fig.tight_layout(); fig.savefig(f'page3_{tag}.png',dpi=90); print('DONE',out,flush=True)
