import numpy as np,json,glob,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1]; th=np.arange(360.)
ostat=dict(rms=float(np.sqrt(np.mean(obs**2))),range=float(obs.max()-obs.min()))
Mm=0.008532635720463873; res={}
def load(M,R):
    a=np.load(f'h_M{M:.5f}_R{R}.npy'); ps=a[:,0]; full=np.full((360,4),np.nan)
    for p,row in zip(ps,a[:,1:]): full[int(p)]=row; full[int(360-p)%360]=row
    return full
for R in [0.05,0.3,0.6,0.9]:
    for M in (Mm,0.0):
        try: d=load(M,R)
        except Exception: continue
        ok=np.isfinite(d[:,0]); lz=np.log(d[ok,0]); dt=-(lz-lz.mean())   # dT/T ~ -dln(1+z) (exact log form)
        # fit peculiar dipole (ring projection) + monopole to Planck over the visible cap
        t=np.radians(th[ok]); X=np.c_[np.ones(ok.sum()),np.cos(t),np.sin(t)]; c,*_=np.linalg.lstsq(X,obs[ok]-dt,rcond=None)
        resid=obs[ok]-dt-X@c
        res[(M>0,R)]=dict(n_valid_deg=int(ok.sum()),cap_deg=[float(th[ok].min()),float(th[ok].max())],rms=float(np.sqrt(np.mean(dt**2))),range=float(np.ptp(dt)),
          one_plus_z=[float(d[ok,0].min()),float(d[ok,0].max())],emit_before_now=[float(np.nanmin(d[:,2])),float(np.nanmax(d[:,2]))],
          vpec_ring_kms=float(np.hypot(c[1],c[2])*299792.458),vpec_ring_angle=float(np.degrees(np.arctan2(c[2],c[1]))%360),
          resid_rms=float(np.sqrt(np.mean(resid**2))),obs_rms_in_cap=float(np.std(obs[ok])),dt=dt,ok=ok,d=d)
out={'planck_full_ring':ostat}
for R in [0.05,0.3,0.6,0.9]:
    if (True,R) in res and (False,R) in res:
        a,b=res[(True,R)],res[(False,R)]
        out[str(R)]=dict(with_mass={k:v for k,v in a.items() if k not in('dt','ok','d')},no_mass={k:v for k,v in b.items() if k not in('dt','ok','d')},
          rms_change_pct_mass_vs_none=100*(a['rms']/b['rms']-1),range_change_pct=100*(a['range']/b['range']-1),
          range_over_obs_with=a['range']/ostat['range'],range_over_obs_none=b['range']/ostat['range'])
        fig,ax=plt.subplots(2,1,figsize=(13,10),sharex=True,gridspec_kw=dict(height_ratios=[2,1]))
        for lab,r,col in (('with core mass',a,'C3'),('no core mass',b,'C0')):
            x=th[r['ok']]; ax[0].plot(x,r['dt'],col,lw=2,label=f"{lab}: RMS {r['rms']:.4f}, range {r['range']:.4f}, range/obs {r['range']/ostat['range']:.1f}")
            ax[1].plot(x,r['d'][r['ok'],2],col,lw=2,label=lab)
        ax[0].plot(th,obs,'k',lw=1,label=f"Planck ring (full): RMS {ostat['rms']:.5f}, range {ostat['range']:.4f}")
        ax[0].set_ylabel('dT/T (model mean removed over visible directions)');ax[1].set_ylabel('emission time before now (units D/c)');ax[1].set_xlabel('ring angle (deg), CoE at 180')
        ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=1,fontsize=10);ax[1].legend(fontsize=9)
        ax[0].set_title(f'recombination sphere radius {R} D; blank = line of sight never meets plasma',y=1.22,fontsize=11)
        for x in ax: x.grid(alpha=.3); x.set_xlim(0,360)
        fig.tight_layout();fig.savefig(f'v5_Rrec{R}.png',dpi=110);plt.close(fig)
json.dump(out,open('v5_stats.json','w'),indent=1);print(json.dumps(out,indent=1))
