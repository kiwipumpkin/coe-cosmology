import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1];PR=float(np.ptp(obs))
psi=np.arange(0,181,1.0);p=np.radians(psi);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(0.5))/4;w=w/w.sum()
P=json.load(open('../../three_vel/point/prof_c4.json'));S={}
for k in ('4','8','15','30'):
    z=np.load(f'sweep/s{k}_mass.npz')['o'][:,1];assert np.isfinite(z).all();T=1/z;d=T/np.sum(w*T)-1
    X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],d*sw,rcond=None)[0];res=d-X@c
    S[k]=dict(f=P[k]['x'][0],d=d,res=res,rng=np.ptp(d),rrng=np.ptp(res),dip=abs(c[1]),rms=np.sqrt(np.sum(w*d*d)),rrms=np.sqrt(np.sum(w*res*res)))
    s=S[k];print(f"{k}c mass {s['f']:.6f} range {s['rng']:.5f} ({s['rng']/PR:.2f}x) dipole {s['dip']:.5f} RMS {s['rms']:.5f} resid range {s['rrng']:.5f} ({s['rrng']/PR:.2f}x) resid RMS {s['rrms']:.5f}")
A=np.r_[psi,360-psi[-2::-1]];m=lambda y:np.r_[y,y[-2::-1]]
for key,out,yl,tt in (('d','dTT_360_mass_speeds.png','dT/T (relative to full-sky mean)','with mass; green 4.2 Gpc at every speed (D scaled); recombination radius 0.5 D; T ∝ 1/(1+z)'),
                      ('res','dTT_360_mass_speeds_residual.png','dT/T left after removing the dipole','with mass, dipole (cosine fit) removed; green 4.2 Gpc; recombination radius 0.5 D; T ∝ 1/(1+z)')):
    fig,ax=plt.subplots(figsize=(13,8))
    for k,col in zip(S,('C0','C2','C1','C3')):
        s=S[k];ax.plot(A,m(s[key]),color=col,lw=2.5,label=f"{k} c, core mass {s['f']:.4f} of critical: range {s['rng']/PR:.1f}× Planck, left after dipole removal {s['rrng']/PR:.1f}× Planck")
    ax.axhline(0,color='k',lw=.6);ax.grid(alpha=.3);ax.set_xlim(0,360);ax.set_xticks(range(0,361,30))
    ax.set_xlabel('angle on the sky from the anti-CoE direction (0°/360° = away from CoE, 180° = toward CoE; mirrored)');ax.set_ylabel(yl)
    ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10.5,frameon=False);ax.set_title(tt,fontsize=10,pad=-14)
    fig.savefig(out,dpi=110,bbox_inches='tight')
