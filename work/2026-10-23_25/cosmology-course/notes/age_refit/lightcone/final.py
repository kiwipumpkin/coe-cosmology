"""usage: final.py <law>  -> res_<law>.json, page_<law>.png   (6-deg grid only)"""
import sys,json,pickle,os,numpy as np; from lc import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
law=sys.argv[1]; st=json.load(open(f'ck_{law}.json')); x=np.array(st['nm']['x']); fin=st['final6']; l0=fin['l0']; af=np.array(fin['af'])
L,tn=make(law,x); o=obs(L,tn); ells=np.arange(0,360,6.)
rk=f'rays_{law}.pkl'
if os.path.exists(rk): RD=pickle.load(open(rk,'rb'))
else:
    _,RD=cmb(L,tn,ells,tmode='A',budget=250000); pickle.dump(RD,open(rk,'wb'))
lam=o['Hbar']*F.C/(70/F.KMSMPC)            # time unit -> Gyr after rescaling to Hbar=70 (scale symmetry)
G2=lambda t: t/F.C*lam                 # Gpc/c -> Gyr  (1 Gpc/c = 1/0.3066 Gyr), rescaled
out=dict(law=law,names=NAMES[law],x=x.tolist(),chi2=fin['chi2'],sigma=fin['sigma'],l0=l0,nfree=len(x)+1,
  U_c=o['U'],Hpar=70*o['Hpar']/o['Hbar'],Hper=70*o['Hper']/o['Hbar'],H_bar=70.,D_Gpc=o['R']*lam,age_Gyr=G2(tn),
  M_kg=L.P['M']*lam, rs_over_D=float(2*L.m(np.array([L.xmin]))[0]/o['R']),
  m_enc_obs_over_D=float(2*L.m(np.array([0.]))[0]/o['R']))
stops={}
for r in RD:
    d=dict(stop=r['stop'],t_end_over_tn=float(r['t'][-1]/tn) if len(r['t']) else None)
    if len(r['t']): xx=np.array([r['x'][-1]]); Re=L.R(r['t'][-1],xx)[0]; d.update(R_end=float(Re),twoM_over_R_end=float(2*L.m(xx)[0]/Re),ncross=int(r['nc'].max()),x_end=float(r['x'][-1]))
    stops[int(round(np.degrees(r['th'])))]=d
out['rays']=stops
for tm in ('A','B'):
    res,_=cmb(L,tn,ells,tmode=tm,RD=RD); d=dict(K=res.get('K'),tstar_Gyr=G2(res['tstar']) if res.get('tstar') else None)
    if res.get('K'):
        sc=score(res['rows'],l0,ells); d.update(sc); d['rows']=[{k:(float(v) if isinstance(v,(np.floating,float)) else v) for k,v in q.items()} for q in res['rows']]
        d['normalisation_K_needed']=res['K']; d['T_emit_over_Tobs_mean']=3000/2.725
        for q in d['rows']:
            if 't_e' in q: q['t_e_Gyr']=G2(q['t_e']); q['R_e_over_D']=float(L.R(q['t_e'],np.array([q['x_e']]))[0]/o['R']); q['U_e_c']=float(L.Rdot(q['t_e'],np.array([q['x_e']]))[0])
        # observer peculiar velocity variants (radial, |beta|<c): Doppler+aberration on the 6-deg ring, then m=1 removed and rescored
        ok=[q for q in res['rows'] if 'T_obs' in q]
        if len(ok)==len(ells):
            Tt=np.array([q['T_obs'] for q in ok]); pv={}
            for b in (-0.5,-0.2,0.2,0.5):
                g=1/np.sqrt(1-b*b); cp=np.cos(np.radians(ells)); c=(cp-b)/(1-b*cp); th=np.degrees(np.arccos(np.clip(c,-1,1))); th=np.where(ells>180,360-th,th)
                Ts=np.exp(np.interp(th,np.append(ells,360),np.log(np.append(Tt,Tt[0]))))/(g*(1-b*cp)); Ts*=2.725/Ts.mean()
                rows2=[dict(q,T_obs=float(v)) for q,v in zip(ok,Ts)]; s2=score(rows2,l0,ells)
                pv[str(b)]={k:s2[k] for k in ('dTT_rms_m1rem','resid_rms','cancel_pct','corr_with_Gm1','T_range')}
            d['peculiar_velocity_variants']=pv
    else: d['K_scan']=res.get('K_scan')
    out['cmb_'+tm]=d
json.dump(out,open(f'res_{law}.json','w'),indent=1,default=float)
# page
A=out['cmb_A']; fig,ax=plt.subplots(2,1,figsize=(19,11))
ell=np.arange(360.); ax[0].errorbar(F.IMIG,F.Y,yerr=0.5*(F.MLO+F.MHI),fmt='o',ms=7,c='k',label='Migkas 71 bins')
ax[0].plot(ell,F.shifted(af,l0),lw=4,c='C3',label=f'light-cone exact GR: χ²={out["chi2"]:.1f}, {out["sigma"]:.1f}σ, l0={l0:.0f}°')
ax[0].set_ylabel('A / A_all',fontsize=22); ax[0].legend(fontsize=18)
if 'gal_l' in A:
    gl=np.array(A['gal_l']); i=np.argsort(gl)
    ax[1].plot(gl[i],np.array(A['dTT'])[i],'o-',lw=3,c='C0',label=f'model dT/T (m=1 removed), rms {A["dTT_rms_m1rem"]:.3g}')
    ax[1].plot(gl[i],np.array(A['Gm1'])[i],'-',lw=3,c='k',label='G−1 (m=1 removed)')
    ax[1].plot(gl[i],np.array(A['resid'])[i],':',lw=3,c='C2',label=f'residual {A["resid_rms"]:.3g}, cancel {A["cancel_pct"]:.0f}%, corr {A["corr_with_Gm1"]:+.2f}')
    ns=A['n_total']-A['n_emit']
    if ns: ax[1].text(0.01,0.04,f'{ns} of {A["n_total"]} directions flagged (no emitter / integrator stopped), not scored',transform=ax[1].transAxes,fontsize=16)
ax[1].set_xlabel('galactic longitude ℓ (deg)',fontsize=22); ax[1].set_ylabel('dT/T',fontsize=22); ax[1].legend(fontsize=16)
for a in ax: a.set_xlim(0,360); a.tick_params(labelsize=18); a.grid()
fig.suptitle(f'{law}: U/c={out["U_c"]:.2f}  D={out["D_Gpc"]:.2f} Gpc  M={out["M_kg"]:.2g} kg  age={out["age_Gyr"]:.1f} Gyr  K needed={A.get("K") or float("nan"):.3g} K  T {A.get("T_range",["?","?"])}',fontsize=18)
fig.tight_layout(); fig.savefig(f'page_{law}.png',dpi=90); print('DONE',law,flush=True)
