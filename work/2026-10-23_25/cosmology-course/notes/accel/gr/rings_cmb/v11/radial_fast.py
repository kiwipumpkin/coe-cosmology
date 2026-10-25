import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
import sys;u=float(sys.argv[1]);Rrec=0.5;t0=1/u;tc=Rrec/u;G=4.2/(t0-tc);Gy=G*3.26156
a=np.radians(np.arange(0,360,2.0));xu=u*tc;rho=0.4*Rrec
xc,yc=xu+rho*np.cos(a),rho*np.sin(a);r=np.hypot(xc,yc);v=r/tc;xt,yt=xc*t0/tc,yc*t0/tc
print(f'D={G:.4f} Gpc, today {t0*Gy:.3f} Gyr, crossing {tc*Gy:.3f} Gyr, ring radius {rho*G:.3f} Gpc')
for d in (0,90,180):
    i=d//2;print(f'{d} deg: dist from CoE at crossing {r[i]*G:.3f} Gpc, speed {v[i]:.4f} c, distance travelled from origin by today {v[i]*t0*G:.3f} Gpc ({v[i]*(t0-tc)*G:.3f} Gpc since crossing)')
X=np.c_[2*xt,2*yt,np.ones_like(xt)];c=np.linalg.lstsq(X,xt**2+yt**2,rcond=None)[0];R=np.sqrt(c[2]+c[0]**2+c[1]**2);dev=np.abs(np.hypot(xt-c[0],yt-c[1])-R).max()
print(f'red today circle fit: centre ({c[0]*G:.4f},{c[1]*G:.1e}) Gpc (us today at {u*t0*G:.4f}), radius {R*G:.4f} Gpc, max dev {dev*G:.1e} Gpc')
dmin=np.hypot(xt-u*t0,yt).min();print(f'closest red to CoE today {np.hypot(xt,yt).min()*G:.3f} Gpc vs Rrec {Rrec*G:.3f}; our speed {u} c');print(f'green {(t0-tc)*G:.4f} Gpc, min red-to-us today {dmin*G:.4f} -> inside {dmin>=t0-tc}; red radius growth {(R-rho)/(t0-tc):.4f} c')
fig,ax=plt.subplots(figsize=(13,11));th=np.linspace(0,2*np.pi,721)
for i in range(0,180,10): ax.plot([0,xt[i]*G],[0,yt[i]*G],color='gray',lw=.6);ax.plot([xt[i]*G,u*t0*G],[yt[i]*G,0],color='darkorange',lw=1)
for d in (0,90,180): i=d//2;print(f'ray {d}: length {np.hypot(xt[i]-u*t0,yt[i])*G:.4f} Gpc')
ax.plot(xc*G,yc*G,'o',color='lightcoral',ms=4,label=f'particles at our crossing ({tc*Gy:.1f} Gyr): ring centred on us, radius 0.4 × recombination radius = {rho*G:.2f} Gpc')
ax.plot(xt*G,yt*G,'o',color='darkred',ms=4,label=f'same particles today ({t0*Gy:.1f} Gyr), each flown straight out from the CoE at its own constant speed (gray lines); orange = straight light rays to us')
ax.plot((u*t0+(t0-tc)*np.cos(th))*G,(t0-tc)*np.sin(th)*G,'g-',lw=2,label=f'how far we can see: c × time since we crossed = {(t0-tc)*G:.2f} Gpc, centred on us today')
ax.fill(Rrec*G*np.cos(th),Rrec*G*np.sin(th),color='lightblue',alpha=.35);ax.plot(Rrec*G*np.cos(th),Rrec*G*np.sin(th),'b-',lw=1.5,label=f'recombination sphere 0.5 D = {Rrec*G:.2f} Gpc (shaded = still plasma)')
for d,off in ((0,(5,8)),(90,(5,8)),(180,(-10,-30))):
    i=d//2;ax.annotate(f'{d}°: {v[i]:.3f} c',(xt[i]*G,yt[i]*G),off,textcoords='offset points',fontsize=10,color='darkred')
ax.plot(xu*G,0,'ko',mfc='none',ms=9);ax.annotate('us at crossing',(xu*G,0),(6,-14),textcoords='offset points')
ax.plot(u*t0*G,0,'k*',ms=15);ax.annotate('us today',(u*t0*G,0),(6,8),textcoords='offset points')
ax.plot(0,0,'k+',ms=18,mew=2);ax.annotate('CoE (origin)',(0,0),(-30,10),textcoords='offset points')
ax.set_aspect('equal');ax.grid(alpha=.3);ax.set_xlabel('Gpc (anti-CoE side to the right)');ax.set_ylabel('Gpc')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
ax.set_title(f'No mass; our speed {u:g} c; D = {G:.2f} Gpc; speeds in units of c',fontsize=11,pad=-14)
fig.savefig(f'red_radial_fast_u{u:g}.png',dpi=110,bbox_inches='tight')
