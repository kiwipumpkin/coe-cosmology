import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1]; th=np.arange(360.); C=299792.458
ost=dict(rms=float(np.sqrt(np.mean(obs**2))),range=float(np.ptp(obs)))
def load(M,fr,P):
    a=np.r_[np.load(f'p_M{M}_f{fr}_p{P}_0.npy'),np.load(f'p_M{M}_f{fr}_p{P}_90.npy')];F=np.full((360,4),np.nan)
    for r in a: F[int(r[0])%360]=r[1:];F[int(360-r[0])%360]=r[1:]
    return F
out={'planck_full_ring':ost}
for fr,P in (('0.5','0.05'),('0.5','0.2'),('0.5','0.5'),('0.2','0.2')):
    R={}
    for lab,M in (('with',"0.00853"),('none',"0.00000")):
        d=load(M,fr,P);ok=np.isfinite(d[:,0]);lz=np.log(d[ok,0]);dt=-(lz-lz.mean())
        t=np.radians(th);X=np.c_[np.ones(360),np.cos(t),np.sin(t)];mod=np.interp(th,th[ok],dt)
        c,*_=np.linalg.lstsq(X,obs-mod,rcond=None);res=obs-mod-X@c
        R[lab]=dict(valid=int(ok.sum()),rms=float(np.sqrt(np.mean(dt**2))),range=float(np.ptp(dt)),one_plus_z=[float(d[ok,0].min()),float(d[ok,0].max())],
          T_obs_K=[3000/float(d[ok,0].max()),3000/float(d[ok,0].min())],emit_proper_time=[float(np.nanmin(d[:,1])),float(np.nanmax(d[:,1]))],
          emit_before_now=[float(np.nanmin(d[:,2])),float(np.nanmax(d[:,2]))],vpec_ring_kms=float(np.hypot(c[1],c[2])*C),vpec_angle=float(np.degrees(np.arctan2(c[2],c[1]))%360),
          resid_rms=float(np.sqrt(np.mean(res**2))),resid_range=float(np.ptp(res)),_d=d,_dt=mod)
    w,n=R['with'],R['none'];pct=lambda a,b:100*(a/b-1) if b>1e-12 else None
    out[fr+'_p'+P]=dict(with_mass={k:v for k,v in w.items() if k[0]!='_'},no_mass={k:v for k,v in n.items() if k[0]!='_'},
       rms_pct_mass_vs_none=pct(w['rms'],n['rms']),range_pct=pct(w['range'],n['range']),range_over_planck_with=w['range']/ost['range'],range_over_planck_none=n['range']/ost['range'])
    fig,ax=plt.subplots(2,1,figsize=(13,10),sharex=True,gridspec_kw=dict(height_ratios=[2,1]))
    for lab,r,col in (('with core mass',w,'C3'),('no core mass',n,'C0')):
        ax[0].plot(th,r['_dt'],col,lw=2,label=f"{lab}: RMS {r['rms']:.4g}, range {r['range']:.4g}, range/Planck {r['range']/ost['range']:.3g}")
        ax[1].plot(th,r['_d'][:,2],col,lw=2,label=lab)
    ax[0].plot(th,obs,'k',lw=1,label=f"Planck ring (full): RMS {ost['rms']:.5f}, range {ost['range']:.4f}")
    ax[0].set_title(f'recombination at {fr} of no-mass age; initial temperature gradient power {P} (hotter toward CoE)',y=1.2);ax[0].set_ylabel('dT/T (mean removed)');ax[1].set_ylabel('emission time before now (our proper time, units D/c)');ax[1].set_xlabel('ring angle (deg), CoE at 180')
    ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=10);ax[1].legend(fontsize=9)
    for x in ax:x.grid(alpha=.3);x.set_xlim(0,360)
    fig.tight_layout();fig.savefig(f'v7_rec{fr}_grad{P}.png',dpi=110);plt.close(fig)
json.dump(out,open('v7_stats.json','w'),indent=1);print(json.dumps(out,indent=1))
