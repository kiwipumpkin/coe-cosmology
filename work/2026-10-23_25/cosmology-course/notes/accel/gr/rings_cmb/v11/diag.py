import numpy as np,v11,base as V,matplotlib,json;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import brentq
M=0.008532635720463873;R=0.5;F=v11.field(M,R);wo,t0,T,W,S=F
tx=brentq(lambda t:V.Rfun(t,wo,M)-R,3e-3,t0)
X=lambda w:V.Rfun(t0,w,M)        # comoving label drawn at the distance that matter has today
times=[tx+0.25*(t0-tx)*k for k in (1,2,3)]+[t0];bend={}
fig,ax=plt.subplots(2,2,figsize=(14,14));ax=ax.ravel()
for i,to in enumerate(times):
    o,P=v11.trace(F,R,np.arange(0,360,5.0),to,keep=True);a=ax[i]
    wr=brentq(lambda w:V.Rfun(to,w,M)-R,W[0]+1e-6,40);a.add_patch(plt.Circle((0,0),X(wr),color='orange',alpha=.35,label=f'plasma at t={to:.3f} (matter not yet crossed {R} D)'))
    a.plot(0,0,'k*',ms=14,label='CoE');a.plot(X(wo),0,'bo',ms=8,label=f'us (fixed comoving position), t={to:.3f} D/c')
    a.plot([],[],'bx',label=f'we crossed recombination at t_x={tx:.3f}')
    ex=[X(w)*np.cos(p) for w,p in zip(o[:,4],o[:,5])];ey=[X(w)*np.sin(p) for w,p in zip(o[:,4],o[:,5])]
    a.plot(np.r_[ex,ex[:1]],np.r_[ey,ey[:1]],'r-',lw=2,label='big bubble edge: last-scattering points')
    sx,sy=[],[];dev=0
    for (ps,tt,ww,pp) in P:
        xs=np.array([X(w) for w in ww])*np.cos(pp);ys=np.array([X(w) for w in ww])*np.sin(pp)
        if int(ps)%30==0 and True: a.plot(xs,ys,'k-',lw=.6)
        if tt[-1]<tx: j=np.argmin(abs(tt-tx));sx.append(xs[j]);sy.append(ys[j])
        else: sx.append(xs[-1]);sy.append(ys[-1])   # ray entered plasma after t_x: clip at the plasma boundary (its last-scattering point)
        # bending: max distance of comoving-label path (w-plane) from straight chord
        qx,qy=ww*np.cos(pp),ww*np.sin(pp);d=np.abs((qx[-1]-qx[0])*(qy-qy[0])-(qy[-1]-qy[0])*(qx-qx[0]))/max(np.hypot(qx[-1]-qx[0],qy[-1]-qy[0]),1e-12);dev=max(dev,d.max())
    if sx: a.plot(np.r_[sx,sx[:1]],np.r_[sy,sy[:1]],'g--',lw=2,label='small bubble: points whose light, sent at t_x, reaches us now')
    a.plot([],[],'k-',lw=.6,label='sample light rays (comoving)');bend[f'{to:.3f}']=float(dev)
    a.set_aspect('equal');a.grid(alpha=.3);a.set_title(f'observation time {to:.3f} D/c'+(' (today)' if to==t0 else ''));a.legend(fontsize=7,loc='lower left')
fig.suptitle(f'v11 comoving ray diagram, with core mass, R_rec={R} D. Each point drawn at the distance its matter has TODAY (units of D)')
fig.tight_layout();fig.savefig('v11_ray_diagram.png',dpi=100)
o,_=v11.trace(F,R,[0,30,150,180]);json.dump(dict(tx=tx,max_bending_wplane=bend,emit_points_today_dist=[[float(p),float(X(w)),float(np.degrees(ph))] for p,w,ph in zip(o[:,0],o[:,4],o[:,5])]),open('diag.json','w'),indent=1)
print(open('diag.json').read())
