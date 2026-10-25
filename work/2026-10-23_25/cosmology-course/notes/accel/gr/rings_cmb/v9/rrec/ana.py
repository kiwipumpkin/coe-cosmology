import numpy as np,json,glob,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../../three_vel/point/cmb/obs_ring.npy')[:-1];th=np.arange(360.);C=299792.458
ost=dict(rms=float(np.sqrt(np.mean(obs**2))),range=float(np.ptp(obs)))
def load(M,R):
    fs=sorted(glob.glob(f'h_M{M}_R{R}_*.npy')) or [f'../../v5/h_M{M}_R{R}.npy'];a=np.concatenate([np.load(f) for f in fs]);a=a[np.argsort(a[:,0])]
    ok=np.isfinite(a[:,1]);g=a[ok]
    if len(g)<3: return None
    lo,hi=g[:,0].min(),180.;cov=(th>=lo)&(th<=360-lo)
    ps=np.r_[g[:,0],360-g[::-1,0]];F=lambda j:np.interp(th,ps,np.r_[g[:,j],g[::-1,j]])
    oz=F(1);lz=np.log(oz[cov]);dt=np.full(360,np.nan);dt[cov]=-(lz-lz.mean())
    t=np.radians(th[cov]);X=np.c_[np.ones(cov.sum()),np.cos(t),np.sin(t)];c,*_=np.linalg.lstsq(X,obs[cov]-dt[cov],rcond=None);res=obs[cov]-dt[cov]-X@c
    te=F(3);te[~cov]=np.nan
    return dict(cov=[float(lo),float(360-lo)],cov_deg=int(cov.sum()),rms=float(np.sqrt(np.nanmean(dt**2))),range=float(np.nanmax(dt)-np.nanmin(dt)),
      one_plus_z=[float(g[:,1].min()),float(g[:,1].max())],emit_before_now=[float(g[:,3].min()),float(g[:,3].max())],
      vpec_kms=float(np.hypot(c[1],c[2])*C),vpec_angle=float(np.degrees(np.arctan2(c[2],c[1]))%360),resid_rms=float(np.sqrt(np.mean(res**2))),
      planck_rms_same_dirs=float(np.std(obs[cov])),_dt=dt,_te=te,_g=g)
out={'planck':ost};RR=[0.05,0.3,0.6,0.9,0.95,0.99];best=None
for R in RR:
    w,n=load('0.00853',R),load('0.00000',R)
    if w is None or n is None: continue
    st=dict(with_mass={k:v for k,v in w.items() if k[0]!='_'},no_mass={k:v for k,v in n.items() if k[0]!='_'},
      rms_pct_mass_vs_none=100*(w['rms']/n['rms']-1),range_pct=100*(w['range']/n['range']-1),range_over_planck_with=w['range']/ost['range'],range_over_planck_none=n['range']/ost['range'])
    out[str(R)]=st
    if best is None or w['resid_rms']<out[str(best)]['with_mass']['resid_rms']: best=R
    fig,ax=plt.subplots(3,1,figsize=(13,13),gridspec_kw=dict(height_ratios=[2,1,1]))
    for lab,r,col in (('with core mass',w,'C3'),('no core mass',n,'C0')):
        ax[0].plot(th,r['_dt'],col,lw=2,label=f"{lab}: RMS {r['rms']:.4g}, range {r['range']:.4g}, range/Planck {r['range']/ost['range']:.3g}, CMB only {r['cov'][0]:.0f}-{r['cov'][1]:.0f} deg")
        ax[1].plot(th,r['_te'],col,lw=2,label=lab)
        g=r['_g'];ax[2].plot(g[:,4],t0:=g[:,2],col+'o',ms=3,label=f'{lab}: matter crossing radius {R} D (crossing time vs launch label)')
    ax[0].plot(th,obs,'k',lw=1,label=f"Planck ring (full): RMS {ost['rms']:.5f}, range {ost['range']:.4f}")
    ax[0].set_title(f'recombination radius {R} D around the CoE; blank = no CMB in that direction',y=1.24)
    ax[0].set_ylabel('dT/T (mean removed over covered directions)');ax[1].set_ylabel('emission before now (D/c)');ax[1].set_xlabel('ring angle (deg), CoE at 180')
    ax[2].set_xlabel('shell speed label w (larger = farther from CoE)');ax[2].set_ylabel('time it crossed the radius (since start, D/c)')
    ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=9);ax[1].legend(fontsize=9);ax[2].legend(fontsize=8)
    for x in ax[:2]:x.set_xlim(0,360)
    for x in ax:x.grid(alpha=.3)
    fig.tight_layout();fig.savefig(f'v9_Rrec{R}.png',dpi=105);plt.close(fig)
out['best_Rrec_by_resid_with_mass']=best;json.dump(out,open('v9_rrec_stats.json','w'),indent=1)
for R in RR:
    if str(R) in out:
        s=out[str(R)];w,n=s['with_mass'],s['no_mass'];print(R,'cov',w['cov'],n['cov'],'rms',round(w['rms'],4),round(n['rms'],4),'rng',round(w['range'],4),round(n['range'],4),'pct',round(s['rms_pct_mass_vs_none'],1),round(s['range_pct'],1),'r/P',round(s['range_over_planck_with']),round(s['range_over_planck_none']),'z',[round(x,3) for x in w['one_plus_z']],'emit',[round(x,3) for x in w['emit_before_now']],[round(x,3) for x in n['emit_before_now']],'vp',round(w['vpec_kms']),round(n['vpec_kms']),'res',round(w['resid_rms'],4),round(n['resid_rms'],4),'Prms',round(w['planck_rms_same_dirs'],5))
print('best',best)
