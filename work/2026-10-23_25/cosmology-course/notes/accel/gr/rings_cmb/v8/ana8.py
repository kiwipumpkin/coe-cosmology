import numpy as np,json,glob,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1];th=np.arange(360.);C=299792.458
ost=dict(rms=float(np.sqrt(np.mean(obs**2))),range=float(np.ptp(obs)))
def full(a):  # a: rows ps,1+z,tau_e,before_now,(w),R  on 0..180 -> 360
    a=a[np.isfinite(a[:,1])];ps=np.r_[a[:,0],360-a[::-1,0]];
    cols=[np.interp(th,ps,np.r_[a[:,j],a[::-1,j]]) for j in range(1,a.shape[1])];return np.array(cols).T,len(a)
def stats(d):
    lz=np.log(d[:,0]);dt=-(lz-lz.mean());t=np.radians(th);X=np.c_[np.ones(360),np.cos(t),np.sin(t)]
    c,*_=np.linalg.lstsq(X,obs-dt,rcond=None);res=obs-dt-X@c
    return dt,dict(rms=float(np.sqrt(np.mean(dt**2))),range=float(np.ptp(dt)),one_plus_z=[float(d[:,0].min()),float(d[:,0].max())],
      emit_before_now=[float(d[:,2].min()),float(d[:,2].max())],emit_R=[float(d[:,-1].min()),float(d[:,-1].max())],
      vpec_kms=float(np.hypot(c[1],c[2])*C),vpec_angle=float(np.degrees(np.arctan2(c[2],c[1]))%360),resid_rms=float(np.sqrt(np.mean(res**2))))
out={'planck':ost};Mm="0.00853"
for fr,q in (('0.5','0.05'),('0.5','0.2'),('0.05','0.1'),('0.01','0.1')):
    R={}
    for lab,M in (('with',Mm),('none','0.00000')):
        a=np.r_[np.load(f'p_M{M}_f{fr}_q{q}_0.npy'),np.load(f'p_M{M}_f{fr}_q{q}_90.npy')]
        nf=int((a[:,4]<0).sum()) if lab=='with' else 0; nmiss=int((~np.isfinite(a[:,1])).sum())
        d,n=full(a[:,[0,1,2,3,4,5]]);dt,s=stats(d);s.update(directions_traced=len(a),failed=nmiss,hit_fallback_matter=nf);R[lab]=(d,dt,s)
    w,n=R['with'][2],R['none'][2]
    out[f'rec{fr}_q{q}']=dict(with_mass=w,no_mass=n,rms_pct=100*(w['rms']/n['rms']-1),range_pct=100*(w['range']/n['range']-1),
        range_over_planck_with=w['range']/ost['range'],range_over_planck_none=n['range']/ost['range'])
    fig,ax=plt.subplots(3,1,figsize=(13,13),gridspec_kw=dict(height_ratios=[2,1,1]))
    for lab,col,nm in (('with','C3','with core mass'),('none','C0','no core mass')):
        d,dt,s=R[lab];ax[0].plot(th,dt,col,lw=2,label=f"{nm}: RMS {s['rms']:.4g}, range {s['range']:.4g}, range/Planck {s['range']/ost['range']:.3g}")
        ax[1].plot(th,d[:,2],col,lw=2,label=nm)
    ax[0].plot(th,obs,'k',lw=1,label=f"Planck ring (full): RMS {ost['rms']:.5f}, range {ost['range']:.4f}")
    ts=float(fr)*0.950762655447463;r=np.linspace(0.05,3,300);ax[2].plot(r,ts*r**(-float(q)),'k',lw=2,label='emission time (since start) vs CoE distance, same with and without mass')
    ax[0].set_title(f'emission time set by CoE distance only: t_e = {ts:.4g} x (r/D)^-{q}',y=1.2)
    ax[0].set_ylabel('dT/T (mean removed)');ax[1].set_ylabel('emission before now (D/c)');ax[1].set_xlabel('ring angle (deg), CoE at 180')
    ax[2].set_xlabel('distance from CoE (units of D)');ax[2].set_ylabel('emission time since start (D/c)')
    ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=10);ax[1].legend(fontsize=9);ax[2].legend(fontsize=9)
    for x in ax[:2]:x.set_xlim(0,360)
    for x in ax:x.grid(alpha=.3)
    fig.tight_layout();fig.savefig(f'v8_rec{fr}_q{q}.png',dpi=105);plt.close(fig)
# exact no-mass at realistic redshift
for f in sorted(glob.glob('milne_f*.npy')):
    a=np.load(f);d,_=full(a);dt,s=stats(d);out['exact_nomass_'+f[6:-4]]=dict(s,range_over_planck=s['range']/ost['range'])
json.dump(out,open('v8_stats.json','w'),indent=1)
for k,v in out.items():
    if k=='planck':continue
    if 'with_mass' in v:
        w,n=v['with_mass'],v['no_mass'];print(k,'rms',round(w['rms'],5),round(n['rms'],5),'rng',round(w['range'],4),round(n['range'],4),'pct',round(v['rms_pct'],1),round(v['range_pct'],1),'r/P',round(v['range_over_planck_with'],1),round(v['range_over_planck_none'],1),'z',[round(x,2) for x in w['one_plus_z']],[round(x,2) for x in n['one_plus_z']],'emit',[round(x,3) for x in w['emit_before_now']],[round(x,3) for x in n['emit_before_now']],'vp',round(w['vpec_kms']),round(w['vpec_angle']),round(n['vpec_kms']),round(n['vpec_angle']),'res',round(w['resid_rms'],5),round(n['resid_rms'],5),'fail',w['failed'],n['failed'],'fallback',w['hit_fallback_matter'],w['directions_traced'])
    else: print(k,{a:(round(b,5) if isinstance(b,float) else b) for a,b in v.items()})
