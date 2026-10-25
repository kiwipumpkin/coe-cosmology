import json,glob,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
rows_out=[]
for f in sorted(glob.glob('gr2/*_*_*.json')):
    d=json.load(open(f)); R=[r for r in d['rows'] if r.get('z_fwd') is not None]
    if not R: rows_out.append((d['tag'],'no direction reaches recombined matter')); continue
    el=(np.array([r['ell'] for r in R])+d['l0'])%360; i=np.argsort(el); el=el[i]
    S=1/(1+np.array([r['z_fwd'] for r in R]))[i]; dtt=S/S.mean()-1
    th=np.radians(el); X=np.c_[np.ones(len(th)),np.cos(th),np.sin(th)]; c=np.linalg.lstsq(X,np.log(S),rcond=None)[0]; Sn=S/np.exp(X[:,1:]@c[1:]); dn=Sn/Sn.mean()-1
    miss=[(np.array(r['ell'])+d['l0'])%360 for r in d['rows'] if r.get('z_fwd') is None]
    fig,ax=plt.subplots(2,1,figsize=(16,11))
    ax[0].plot(el,dtt*1e2,'o-',lw=3,label='dT/T raw (×100)'); ax[0].plot(el,dn*1e2,'s-',lw=3,label='dT/T, m=1 removed (×100, fit on visible part)')
    for m in miss: ax[0].axvspan(m-3,m+3,color='0.8')
    ax[0].set_ylabel('dT/T ×100',fontsize=18); ax[0].legend(fontsize=14); ax[0].set_title(f"{d['tag']}: grey = no recombined matter on our past cone (ray enters trapped region 2Gm/(c²R)>1+2E)",fontsize=15)
    Re=np.array([r['R_e'] for r in R])[i]; ax[1].plot(el,Re/ (d['D_Gpc']/d['D_Gpc']) ,'o-',lw=3,label='emitter radius R_e (code units, D_code≈1)')
    ax[1].set_ylabel('R_e at t_rec',fontsize=18); ax[1].set_xlabel('galactic longitude ℓ (deg)',fontsize=18); ax[1].legend(fontsize=14)
    ax[1].set_title(f"emission time t_rec common; flight time {d['flight_time_Gyr']:.3g} Gyr all directions; <z>={np.mean(1/S-1):.4g}",fontsize=15)
    for a in ax: a.set_xlim(0,360); a.tick_params(labelsize=14); a.grid()
    fig.tight_layout(); fig.savefig(f"gr2/cmbpartial_{d['tag']}.png",dpi=80); plt.close(fig)
    rows_out.append((d['tag'],f"visible {len(R)}/{d['n_dir']}, dT/T raw p2p {np.ptp(dtt):.3g}, m1-removed p2p {np.ptp(dn):.3g}, z range {1/S.max()-1:.4g}-{1/S.min()-1:.4g}"))
for r in rows_out: print(*r)
