import numpy as np,json,pickle,glob,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1];th=np.arange(360.);C=299792.458
ost=dict(rms=float(np.sqrt(np.mean(obs**2))),range=float(np.ptp(obs)))
D={}
for M in ('0.00853','0.00000'):
    for f in glob.glob(f'v10_M{M}_*.pkl'):
        d=pickle.load(open(f,'rb'));D.setdefault(M,dict(t0=d['t0'],wo=d['wo'],res={}))['res'].update(d['res'])
t=np.radians(th);X=np.c_[np.ones(360),np.cos(t),np.sin(t)]
def one(o,t0,wo,R):
    full=np.zeros((360,o.shape[1]))
    for r in o: full[int(r[0])%360]=r;full[int(360-r[0])%360]=r
    dt=-np.log(full[:,1]);dt-=dt.mean();c,*_=np.linalg.lstsq(X,obs-dt,rcond=None);res=obs-dt-X@c
    tx=float(np.interp(R,[np.exp(0)]*0+[0],[0])) if False else None
    return full,dt,dict(coverage_deg=int(np.isfinite(full[:,1]).sum()),rms=float(np.sqrt(np.mean(dt**2))),range=float(np.ptp(dt)),one_plus_z=[float(full[:,1].min()),float(full[:,1].max())],
        emit_time=[float(full[:,2].min()),float(full[:,2].max())],emit_before_now=[float(full[:,3].min()),float(full[:,3].max())],
        vpec_kms=float(np.hypot(c[1],c[2])*C),vpec_angle=float(np.degrees(np.arctan2(c[2],c[1]))%360),resid_rms=float(np.sqrt(np.mean(res**2))))
out={'planck':ost};Rs=sorted(set(D['0.00853']['res'])&set(D['0.00000']['res']))
for R in Rs:
    A={}
    for M in D: A[M]=one(D[M]['res'][R][0],D[M]['t0'],D[M]['wo'],R)
    w,n=A['0.00853'][2],A['0.00000'][2]
    out[str(R)]=dict(with_mass=w,no_mass=n,rms_pct=100*(w['rms']/n['rms']-1),range_pct=100*(w['range']/n['range']-1),range_over_planck_with=w['range']/ost['range'],range_over_planck_none=n['range']/ost['range'])
    fig,ax=plt.subplots(2,1,figsize=(13,10),sharex=True,gridspec_kw=dict(height_ratios=[2,1]))
    for M,col,lab in (('0.00853','C3','with core mass'),('0.00000','C0','no core mass')):
        f,dt,s=A[M];ax[0].plot(th,dt,col,lw=2,label=f"{lab}: RMS {s['rms']:.4g}, range {s['range']:.4g}, range/Planck {s['range']/ost['range']:.3g}, coverage {s['coverage_deg']}/360 deg")
        ax[1].plot(th,f[:,2],col,lw=2,label=lab)
    ax[0].plot(th,obs,'k',lw=1,label=f"Planck ring (full): RMS {ost['rms']:.5f}, range {ost['range']:.4f}")
    ax[0].set_title(f"Chris kinematics; recombination radius {R} D",y=1.2);ax[0].set_ylabel('dT/T (mean removed)')
    ax[1].set_ylabel('emission time since start (D/c)');ax[1].set_xlabel('ring angle (deg), CoE at 180')
    ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=9);ax[1].legend(fontsize=9)
    for x in ax:x.grid(alpha=.3);x.set_xlim(0,360)
    fig.tight_layout();fig.savefig(f'v10_Rrec{R}.png',dpi=105);plt.close(fig)
best=min(Rs,key=lambda R:out[str(R)]['with_mass']['resid_rms']);out['best_Rrec_with_mass']=best
out['best_Rrec_no_mass']=min(Rs,key=lambda R:out[str(R)]['no_mass']['resid_rms'])
json.dump(out,open('v10_stats.json','w'),indent=1)
for R in Rs:
    s=out[str(R)];w,n=s['with_mass'],s['no_mass'];print(R,'cov',w['coverage_deg'],n['coverage_deg'],'rms',round(w['rms'],4),round(n['rms'],4),'rng',round(w['range'],4),round(n['range'],4),'pct',round(s['rms_pct'],1),round(s['range_pct'],1),'r/P',round(s['range_over_planck_with']),round(s['range_over_planck_none']),'z',[round(x,3) for x in w['one_plus_z']],[round(x,3) for x in n['one_plus_z']],'te',[round(x,4) for x in w['emit_time']],[round(x,4) for x in n['emit_time']],'vp',round(w['vpec_kms']),round(n['vpec_kms']),'res',round(w['resid_rms'],4),round(n['resid_rms'],4))
print('best',best,out['best_Rrec_no_mass'])
