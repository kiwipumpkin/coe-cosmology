import numpy as np,v10,v5 as V,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import brentq
M=0.008532635720463873;R=0.5
F=v10.field(M);wo,t0,U,T,S=F
tx=brentq(lambda t:V.Rfun(t,wo,M)-R,1e-4,t0)
times=[tx+0.25*(t0-tx)*k for k in (1,2,3)]+[t0]
fig,ax=plt.subplots(2,2,figsize=(14,14));ax=ax.ravel()
for i,to in enumerate(times):
    o,p,_,_=v10.rays(M,R,np.arange(0,361,10.0),(wo,to,U,T,S))
    a=ax[i];ro=V.Rfun(to,wo,M)
    a.add_patch(plt.Circle((0,0),R,color='orange',alpha=.35,label=f'plasma (inside recombination radius {R} D)'))
    a.plot(0,0,'k*',ms=14,label='CoE');a.plot(ro,0,'bo',ms=8,label=f'us at t={to:.3f} D/c (r={ro:.3f} D)')
    ts=np.linspace(0,to,100);a.plot([V.Rfun(t,wo,M) for t in ts],0*ts,'b-',lw=1,label='our path')
    a.plot(R,0,'bx',ms=12,mew=3,label=f'our crossing, t_x={tx:.3f}')
    ph=np.linspace(0,2*np.pi,200);rr=to-tx;a.plot(ro+rr*np.cos(ph),rr*np.sin(ph),'g--',lw=2,label=f'small bubble: c(t - t_x) = {rr:.3f} D')
    xs=o[:,5];re=np.array([V.Rfun(te,w,M) for te,w in zip(o[:,2],o[:,4])])
    ex,ey=re*np.cos(xs),re*np.sin(xs);a.plot(np.r_[ex,ex[:1]],np.r_[ey,ey[:1]],'r-',lw=2,label='big bubble edge: last-scattering points (where each ray left plasma)')
    for (ps,tt,rp,thp) in p[::3]: a.plot(rp*np.cos(thp),rp*np.sin(thp),'k-',lw=.6,alpha=.7)
    a.plot([],[],'k-',lw=.6,label='sample light rays (traced back from us)')
    L=1.6*max(ro,R);a.set_xlim(-L*0.6,L);a.set_ylim(-L*0.8,L*0.8);a.set_aspect('equal');a.grid(alpha=.3)
    a.set_title(f'observation time {to:.3f} D/c'+(' (today)' if to==t0 else ''));a.legend(fontsize=7,loc='upper left')
fig.suptitle('Ray diagram, with core mass, recombination radius 0.5 D (positions in units of D; CoE at origin; drawn at emission positions)')
fig.tight_layout();fig.savefig('v10_ray_diagram.png',dpi=100)
