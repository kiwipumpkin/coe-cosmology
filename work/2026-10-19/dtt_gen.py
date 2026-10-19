import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1];PR=float(np.ptp(obs));print('Planck ring range',PR,'rms',np.sqrt(np.mean(obs**2)))
psi=np.arange(0,181,1.0);p=np.radians(psi);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(0.5))/4  # endpoint caps
w=w/w.sum();S={}
import sys
PRE,LM,OUT,TT=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
for tag,lab in (('nomass','no mass'),('mass',LM)):
    o=np.load(f'{PRE}_{tag}.npz')['o'];z=o[:,1];print(lab,'unresolved',np.arange(181)[~np.isfinite(z)].tolist())
    T=1/z;d=T/np.sum(w*T)-1
    X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],d*sw,rcond=None)[0];res=d-X@c
    S[tag]=dict(lab=lab,d=d,rms=float(np.sqrt(np.sum(w*d*d))),rng=float(np.ptp(d)),dip=float(abs(c[1])),rrms=float(np.sqrt(np.sum(w*res*res))),rrng=float(np.ptp(res)),zr=(z.min(),z.max()))
    s=S[tag];print(f"{lab}: RMS {s['rms']:.5f} range {s['rng']:.5f} dipole {s['dip']:.5f} range/Planck {s['rng']/PR:.1f}; resid RMS {s['rrms']:.5f} range {s['rrng']:.5f} resid range/Planck {s['rrng']/PR:.1f}; 1+z {z.min():.4f}-{z.max():.4f}")
a,b=S['mass'],S['nomass'];print(f"mass vs no mass: RMS {100*(a['rms']/b['rms']-1):+.2f}%, range {100*(a['rng']/b['rng']-1):+.2f}%, resid RMS {100*(a['rrms']/b['rrms']-1):+.2f}%")
fig,ax=plt.subplots(figsize=(13,8))
for tag,col in (('nomass','C0'),('mass','C3')):
    s=S[tag];ax.plot(np.r_[psi,360-psi[-2::-1]],np.r_[s['d'],s['d'][-2::-1]],col,lw=2.5,label=f"{s['lab']}: RMS {s['rms']:.4f}, range {s['rng']:.4f} ({s['rng']/PR:.0f}× Planck ring), dipole {s['dip']:.4f}; after dipole removal RMS {s['rrms']:.4f}, range {s['rrng']:.4f}")
ax.axhline(0,color='k',lw=.6);ax.grid(alpha=.3);ax.set_xlim(0,360);ax.set_xticks(range(0,361,30))
ax.set_xlabel('angle on the sky from the anti-CoE direction (0°/360° = away from CoE, 180° = toward CoE; mirrored)');ax.set_ylabel('dT/T (relative to full-sky mean)')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
ax.set_title(TT,fontsize=10,pad=-14)
fig.savefig(OUT,dpi=110,bbox_inches='tight')
