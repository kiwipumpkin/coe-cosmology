import numpy as np,json,glob,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1];th=np.arange(360.);C=299792.458
ost=dict(rms=float(np.sqrt(np.mean(obs**2))),range=float(np.ptp(obs)))
a=np.concatenate([np.load(f) for f in sorted(glob.glob('p_M0.00853_Z1100_ti0.01_*.npy'))]);a=a[np.argsort(a[:,0])]
fb=a[:,6]==1;fail=a[:,7]==1;good=~fb&~fail&np.isfinite(a[:,1])
g=a[good];ps=np.r_[g[:,0],360-g[::-1,0]];F=lambda j:np.interp(th,ps,np.r_[g[:,j],g[::-1,j]])
oz=F(1);lz=np.log(oz);dt=-(lz-lz.mean());cl=np.cos(np.radians(th))
t=np.radians(th);X=np.c_[np.ones(360),np.cos(t),np.sin(t)]
def fit(m):
    c,*_=np.linalg.lstsq(X,obs-m,rcond=None);r=obs-m-X@c;return float(np.hypot(c[1],c[2])*C),float(np.degrees(np.arctan2(c[2],c[1]))%360),float(np.sqrt(np.mean(r**2)))
vw=fit(dt);vn=fit(np.zeros(360))
t0m,t0n=0.9239972648439513,0.950762655447463;Dw,Dn=4.23/t0m,4.23/t0n
out=dict(planck=ost,with_mass=dict(directions_traced=len(a),clean=int(good.sum()),fallback_hits=a[fb,0].tolist(),fallback_one_plus_z=a[fb,1].tolist(),failed=a[fail,0].tolist(),
  rms=float(np.sqrt(np.mean(dt**2))),range=float(np.ptp(dt)),range_over_planck=float(np.ptp(dt)/ost['range']),one_plus_z=[float(g[:,1].min()),float(g[:,1].max())],
  T_today_K=[3000/float(g[:,1].max()),3000/float(g[:,1].min())],emit_proper_time=[float(g[:,2].min()),float(g[:,2].max())],emit_R=[float(g[:,5].min()),float(g[:,5].max())],
  vpec_kms=vw[0],vpec_angle=vw[1],resid_rms=vw[2],D_Gpc=Dw),
  no_mass=dict(exact='homogeneous Milne: 1+z=1100 every direction, emission at the same proper time 0.000864 for all matter',rms=0.0,range=0.0,vpec_kms=vn[0],vpec_angle=vn[1],resid_rms=vn[2],D_Gpc=Dn,
  numeric_check_psi0=float(np.load('single/p_M0.00000_Z1100_ti0.01_0.npy')[0,1])))
json.dump(out,open('v9_stats.json','w'),indent=1);print(json.dumps(out,indent=1))
fig,ax=plt.subplots(3,1,figsize=(13,13),gridspec_kw=dict(height_ratios=[2,1,1]))
w=out['with_mass']
ax[0].plot(th,dt,'C3',lw=2,label=f"with core mass: RMS {w['rms']:.4g}, range {w['range']:.4g}, range/Planck {w['range_over_planck']:.3g}, 1+z {w['one_plus_z'][0]:.0f}-{w['one_plus_z'][1]:.0f}")
ax[0].plot(th,0*th,'C0',lw=2,label="no core mass (exact): RMS 0, range 0, 1+z = 1100 everywhere")
ax[0].plot(th,obs,'k',lw=1,label=f"Planck ring (full): RMS {ost['rms']:.5f}, range {ost['range']:.4f}")
for x0 in a[fb,0]: 
    for xx in (x0,360-x0): ax[0].axvline(xx,color='0.7',lw=0.8)
ax[0].set_title('computed adiabatic cooling, uniform initial state; grey lines = sight lines ending in fallback matter (excluded)',y=1.22)
ax[1].plot(th,F(3),'C3',lw=2,label='with core mass');ax[1].plot(th,0*th+t0n-0.950762655447463/1100,'C0',lw=2,label='no core mass')
ax[2].plot(g[:,5],g[:,2],'C3o-',label='with core mass (computed emission events)');ax[2].axhline(0.950762655447463/1100,color='C0',lw=2,label='no core mass (all emit at once)')
ax[0].set_ylabel('dT/T (mean removed)');ax[1].set_ylabel('emission before now (D/c)');ax[1].set_xlabel('ring angle (deg), CoE at 180')
ax[2].set_xlabel('distance from CoE at emission (units of D)');ax[2].set_ylabel('emission time since start (D/c)')
ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=10);ax[1].legend(fontsize=9);ax[2].legend(fontsize=9)
for x in ax[:2]:x.set_xlim(0,360)
for x in ax:x.grid(alpha=.3)
fig.tight_layout();fig.savefig('v9_z1100.png',dpi=105)
