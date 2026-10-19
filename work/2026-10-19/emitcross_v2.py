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

fig,axs=plt.subplots(1,2,figsize=(20,9),gridspec_kw=dict(width_ratios=[1,1.1]));th=np.linspace(0,2*np.pi,721)
P=np.r_[phi[ok],-phi[ok][::-1]];RE=np.r_[re[ok],re[ok][::-1]]*G;RT=np.r_[rt[ok],rt[ok][::-1]]*G;X0=xo*G
for i in (0,90,180): print(f' straight ray red->us {i}: {dt[i]*G:.3f} Gpc')
for k,ax in enumerate(axs):
    ax.fill(R*G*np.cos(th),R*G*np.sin(th),color='lightblue',alpha=.3);ax.plot(R*G*np.cos(th),R*G*np.sin(th),'b-',lw=1)
    ax.plot(RE*np.cos(P),RE*np.sin(P),color='mediumblue',lw=3,label=f'blue: emission points on recombination sphere (shaded = plasma), at emission; half-width {np.degrees(np.nanmax(np.abs(phi[ok]))):.1f}° from CoE')
    ax.plot(RT*np.cos(P),RT*np.sin(P),'r-',lw=2.5,label=f'red: same emitters today; straight ray to us 0° {dt[0]*G:.2f}, 90° {dt[90]*G:.2f}, 180° {dt[180]*G:.2f} Gpc')
    ax.plot((xo+rv*np.cos(th))*G,rv*np.sin(th)*G,'g-',lw=2,label=f'green: c × time since we crossed = {rv*G:.2f} Gpc, centred on us today')
    for s_ in (1,-1):
        for i in range(0,181,15):
            if not ok[i]: continue
            c_,sn=np.cos(s_*phi[i]),np.sin(s_*phi[i]);x1,y1,x2,y2=re[i]*G*c_,re[i]*G*sn,rt[i]*G*c_,rt[i]*G*sn
            ax.plot([0,x1],[0,y1],'--',color='gray',lw=.7)
            ax.annotate('',(x2,y2),(x1,y1),arrowprops=dict(arrowstyle='->',color='gray',lw=.9))
            ax.plot([x2,X0],[y2,0],color='darkorange',lw=1.2)
    ax.plot([],[],color='darkorange',lw=1.2,label='orange: straight line from red point to us (every 15°); gray: dashed radial line from CoE, arrow = displacement emission → today')
    ax.plot(X0,0,'k*',ms=15);ax.annotate('us today',(X0,0),(6,8),textcoords='offset points')
    ax.plot(xc*G,0,'ko',mfc='none',ms=8);ax.annotate('us at crossing',(xc*G,0),(6,-14),textcoords='offset points')
    ax.plot(0,0,'k+',ms=18,mew=2);ax.annotate('CoE',(0,0),(-20,10),textcoords='offset points')
    ax.set_aspect('equal');ax.grid(alpha=.3);ax.set_xlabel('Gpc (anti-CoE side to the right)');ax.set_ylabel('Gpc')
X=np.r_[RE*np.cos(P),RT*np.cos(P)];Y=np.r_[RE*np.sin(P),RT*np.sin(P),rv*G,-rv*G];m=1.0
axs[1].set_xlim(X.min()-m,X.max()+m);axs[1].set_ylim(Y.min()-m,Y.max()+m);axs[1].set_title('zoom: arc, us and red',fontsize=11)
axs[0].set_title(f'full view; mass {f:g} of critical; our speed {u:g} c; D = 42 Gpc; today {t0*Gy:.1f} Gyr',fontsize=11)
axs[0].legend(loc='lower left',bbox_to_anchor=(0,1.06),fontsize=10,frameon=False)
fig.savefig(sys.argv[4],dpi=100,bbox_inches='tight')
