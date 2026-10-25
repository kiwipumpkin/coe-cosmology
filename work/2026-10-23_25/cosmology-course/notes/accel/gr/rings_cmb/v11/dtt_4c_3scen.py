import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1];PR=float(np.ptp(obs))
psi=np.arange(0,181,1.0);p=np.radians(psi);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(0.5))/4;w=w/w.sum()
C=[('z4_nomass.npz','4 c, emission at own crossing, no mass (today 34.2 Gyr)','C0'),('z4_m02.npz','4 c, own crossing, 0.2 of critical (today 29.0 Gyr)','C3'),('z4_sim.npz','4 c, simultaneous recombination at 17.12 Gyr, no mass (today 34.2 Gyr)','C2')]
S=[]
for fn,lab,col in C:
    z=np.load(fn)['o'][:,1];assert np.isfinite(z).all();T=1/z;d=T/np.sum(w*T)-1
    X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],d*sw,rcond=None)[0];r=d-X@c
    s=dict(lab=lab,col=col,d=d,r=r,rms=np.sqrt(np.sum(w*d*d)),rng=np.ptp(d),dip=abs(c[1]),rrng=np.ptp(r),rrms=np.sqrt(np.sum(w*r*r)),z=(z.min(),z.max()));S.append(s)
    print(f"{lab}: RMS {s['rms']:.6f} range {s['rng']:.6f} ({s['rng']/PR:.2f}x) dipole {s['dip']:.6f} leftover range {s['rrng']:.3e} ({s['rrng']/PR:.3f}x) leftover RMS {s['rrms']:.3e} 1+z {z.min():.6f}-{z.max():.6f}")
a,b=S[1],S[0];print(f"mass vs no mass: RMS {100*(a['rms']/b['rms']-1):+.2f}%, range {100*(a['rng']/b['rng']-1):+.2f}%")
A=np.r_[psi,360-psi[-2::-1]];m=lambda y:np.r_[y,y[-2::-1]]
for key,out,yl,tt in (('d','dTT_4c_3scenarios.png','dT/T (relative to full-sky mean)','recombination radius 0.5 D, D = 42 Gpc; T ∝ 1/(1+z); equal temperature at recombination'),
                      ('r','dTT_4c_3scenarios_residual.png','dT/T left after removing the dipole (cosine fit)','same three runs, dipole removed')):
    fig,ax=plt.subplots(figsize=(13,8))
    for s in S:
        ax.plot(A,m(s[key]),color=s['col'],lw=2.5,label=(f"{s['lab']}: 1+z = 2.000 in every direction, so dT/T = 0 (range {s['rng']/PR:.1e}× Planck, numerical noise)" if s['rng']<1e-5 else f"{s['lab']}: range {s['rng']/PR:.1f}× Planck, RMS {s['rms']:.4f}, dipole {s['dip']:.4f}, left after dipole removal {s['rrng']/PR:.2f}× Planck"))
    ax.axhline(0,color='k',lw=.6);ax.grid(alpha=.3);ax.set_xlim(0,360);ax.set_xticks(range(0,361,30))
    ax.set_xlabel('angle on the sky from the anti-CoE direction (0°/360° = away from CoE, 180° = toward CoE; mirrored)');ax.set_ylabel(yl)
    ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=9.5,frameon=False);ax.set_title(tt,fontsize=10,pad=-14)
    fig.savefig(out,dpi=110,bbox_inches='tight')
