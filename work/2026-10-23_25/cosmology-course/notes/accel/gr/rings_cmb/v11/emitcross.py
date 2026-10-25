import sys,numpy as np,base as V
u=float(sys.argv[1]);f=float(sys.argv[2]);tag=sys.argv[3];V.uD=u
import v11,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import brentq
M=f*0.5*u**2;R=0.5;F=v11.field(M,R);wo,t0,T,W,S=F
ps=np.arange(0,181,1.0);o,_=v11.trace(F,R,ps);ok=np.isfinite(o[:,1])
te=o[:,2];re=np.array([V.Rfun(a,w,M) for a,w in zip(te,o[:,4])]);rt=np.array([V.Rfun(t0,w,M) for w in o[:,4]]);phi=o[:,5]
xo=V.Rfun(t0,wo,M);tc=brentq(lambda t:V.Rfun(t,wo,M)-R,T[0]*1.01,t0);xc=V.Rfun(tc,wo,M);rv=t0-tc
G=42.0;Gy=G*3.26156
print(f'{tag}: u={u} f={f} M={M:.4f}; found {ok.sum()}/181, missing {ps[~ok].tolist()}; t0={t0:.5f} D/c = {t0*Gy:.3f} Gyr; our crossing {tc*Gy:.3f} Gyr; green {rv*G:.3f} Gpc; our speed today {(V.Rfun(t0*1.0001,wo,M)-V.Rfun(t0*0.9999,wo,M))/(2e-4*t0):.4f} c')
d=lambda r,p:np.hypot(r*np.cos(p)-xo,r*np.sin(p));de,dt=d(re,phi),d(rt,phi)
for i in (0,90,180): print(f' {i}: emit t {te[i]*Gy:.3f} Gyr, light travel {o[i,3]*G:.3f} Gpc, 1+z {o[i,1]:.4f}, red today from us {dt[i]*G:.3f} Gpc, emit pt from us today {de[i]*G:.3f} Gpc, emit radius {re[i]/R:.6f} Rrec')
print(f' delay 0 vs 180: {(te[180]-te[0])*Gy:.3f} Gyr; blue arc half-width from CoE {np.degrees(np.nanmax(np.abs(phi[ok]))):.3f} deg; green inside red: {np.nanmin(dt[ok])>=rv} (closest red {np.nanmin(dt[ok])*G:.3f} Gpc)')
if len(sys.argv)<5: sys.exit()
fig,ax=plt.subplots(figsize=(14,10));th=np.linspace(0,2*np.pi,721)
P=np.r_[phi[ok],-phi[ok][::-1]];RE=np.r_[re[ok],re[ok][::-1]]*G;RT=np.r_[rt[ok],rt[ok][::-1]]*G
ax.fill(R*G*np.cos(th),R*G*np.sin(th),color='lightblue',alpha=.3);ax.plot(R*G*np.cos(th),R*G*np.sin(th),'b-',lw=1)
ax.plot(RE*np.cos(P),RE*np.sin(P),color='mediumblue',lw=3,label=f'blue: emission points on the recombination sphere (21 Gpc, shaded = plasma), at emission; half-width {np.degrees(np.nanmax(np.abs(phi[ok]))):.1f}° from CoE')
ax.plot(RT*np.cos(P),RT*np.sin(P),'r-',lw=2.5,label=f'red: same emitters today, moved radially; from us 0° {dt[0]*G:.1f}, 90° {dt[90]*G:.1f}, 180° {dt[180]*G:.1f} Gpc')
ax.plot((xo+rv*np.cos(th))*G,rv*np.sin(th)*G,'g-',lw=2,label=f'green: c × time since we crossed = {rv*G:.2f} Gpc, centred on us today')
for s in (1,-1):
    for i in range(0,181,15):
        if not ok[i]: continue
        ax.plot([re[i]*G*np.cos(s*phi[i]),rt[i]*G*np.cos(s*phi[i])],[re[i]*G*np.sin(s*phi[i]),rt[i]*G*np.sin(s*phi[i])],color='gray',lw=.6)
for k,(p,tt,ww,pp) in enumerate(v11.trace(F,R,np.arange(0,181,15.0),keep=True)[1]):
    rr=np.array([V.Rfun(a,b,M) for a,b in zip(tt,ww)])
    for s in (1,-1): ax.plot(rr*np.cos(s*pp)*G,rr*np.sin(s*pp)*G,color='darkorange',lw=1,label='orange: light path from emission point to us (every 15°)' if (k==0 and s==1) else None)
ax.plot(xo*G,0,'k*',ms=15);ax.annotate('us today',(xo*G,0),(6,8),textcoords='offset points')
ax.plot(xc*G,0,'ko',mfc='none',ms=8);ax.annotate('us at crossing',(xc*G,0),(6,-14),textcoords='offset points')
ax.plot(0,0,'k+',ms=18,mew=2);ax.annotate('CoE',(0,0),(-20,10),textcoords='offset points')
ax.set_aspect('equal');ax.grid(alpha=.3);ax.set_xlabel('Gpc (anti-CoE side to the right)');ax.set_ylabel('Gpc')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
ax.set_title(f'Mass {f:g} of critical; our speed {u:g} c; D = 42 Gpc; today {t0*Gy:.1f} Gyr',fontsize=11,pad=-14)
fig.savefig(sys.argv[4],dpi=110,bbox_inches='tight')
