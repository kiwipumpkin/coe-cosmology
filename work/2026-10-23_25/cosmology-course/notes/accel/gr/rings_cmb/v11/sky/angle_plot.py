import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1];PR=float(np.ptp(obs))
psi=np.arange(0,181,1.0);p=np.radians(psi);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(0.5))/4;w/=w.sum()
def st(fn):
    z=np.load(fn)['o'][:,1];T=1/z;dd=T/np.sum(w*T)-1
    X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],dd*sw,rcond=None)[0];r=dd-X@c
    return dict(d=dd,r=r,rms=float(np.sqrt(np.sum(w*dd*dd))),rng=float(np.ptp(dd)),dip=float(abs(c[1])),rrng=float(np.ptp(r)))
RUNS=[('distributed/z4_dist_a_own.npz','own crossing, distributed, normalisation a (equal mass)','C0'),
      ('distributed/z4_dist_b_own.npz','own crossing, distributed, normalisation b (Migkas-fitted)','C1'),
      ('distributed/z4_dist_b_sim.npz','simultaneous, distributed, normalisation b (Migkas-fitted)','C2'),
      ('z4_m02.npz','own crossing, point mass, fraction 0.2 (reference)','C3')]
REF=[('LX–T flow',51.1),('Planck hot',82.8),('quadrupole',83.8),('Migkas low-H₀',84.5),('octopole',91.4),('Planck cold',97.2),('YSZ–T flow',99.6),('CatWISE',108.7)]
S=[(st(f),n,c) for f,n,c in RUNS]
a=np.arange(0,361);idx=np.where(a<=180,a,360-a)
for key,out,yl in (('d','dTT_angle_4runs_paper.png','dT/T'),('r','dTT_angle_4runs_paper_residual.png','dT/T left after removing monopole + dipole')):
    fig,ax=plt.subplots(figsize=(16,9))
    for s,n,c in S:
        lab=f"{n}: range {s['rng']/PR:.1f}× Planck, RMS {s['rms']:.4f}, dipole {s['dip']:.4f}, left after dipole {s['rrng']/PR:.2f}× Planck"
        ax.plot(a,s[key][idx],color=c,lw=2.2,label=lab);print(lab)
    yl0,yl1=ax.get_ylim()
    for k,(nm,ang) in enumerate(REF):
        for x in (ang,360-ang):
            ax.axvline(x,color='0.4',ls=':',lw=1)
            ax.text(x,yl1-(0.03+0.11*(k%4))*(yl1-yl0),f' {nm} {ang}°',rotation=90,fontsize=8,va='top',ha='center',color='0.25',bbox=dict(fc='w',ec='none',alpha=.7,pad=0.5))
    ax.set_xlim(0,360);ax.set_xticks(range(0,361,30));ax.axhline(0,color='k',lw=.6)
    ax.set_xlabel('angle from the anti-CoE direction (°)  —  0° = anti-CoE (l, b) = (0°, +15°), 180° = CoE (180°, −15°); paper direction, main.tex lines 359–362')
    ax.set_ylabel(yl);ax.grid(alpha=.25)
    ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.01),fontsize=10.5,frameon=False,title=('dT/T vs angle, speed 4c' if key=='d' else 'dT/T after removing monopole + dipole, speed 4c')+f'; Planck ring range = {PR:.7f}; dotted: reference directions (mirrored at 360° − angle)')
    fig.savefig(out,dpi=90,bbox_inches='tight');plt.close(fig)
