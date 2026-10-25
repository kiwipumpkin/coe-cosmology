import sys,numpy as np,v11,base as V,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import brentq
M=float(sys.argv[1]);R=float(sys.argv[2]);tag=sys.argv[3]
F=v11.field(M,R);wo,t0,T,W,S=F
ps=np.arange(0,181,1.0);o,_=v11.trace(F,R,ps)
ok=np.isfinite(o[:,1]);print('found',ok.sum(),'of',len(ps))
te=o[:,2];re=np.array([V.Rfun(a,w,M) for a,w in zip(te,o[:,4])]);rt=np.array([V.Rfun(t0,w,M) for w in o[:,4]]);phi=o[:,5]
xo=V.Rfun(t0,wo,M);tc=brentq(lambda t:V.Rfun(t,wo,M)-R,T[0]*1.01,t0);xc=V.Rfun(tc,wo,M);rv=t0-tc
d=lambda r,p:np.hypot(r*np.cos(p)-xo,r*np.sin(p))
de,dt=d(re,phi),d(rt,phi)
G=4.58
for i in (0,90,180):
    print(f"sky {ps[i]:.0f}: emit pos x={re[i]*np.cos(phi[i]):.4f} y={re[i]*np.sin(phi[i]):.4f}, emit dist from us {de[i]:.4f} D={de[i]*G:.3f} Gpc, today dist {dt[i]:.4f} D={dt[i]*G:.3f} Gpc, light travel {o[i,3]:.4f} D={o[i,3]*G:.3f} Gpc, t_emit {te[i]:.4f}, 1+z {o[i,1]:.4f}")
print(f"t0 {t0:.4f} us today {xo:.4f}, cross t {tc:.4f} at {xc:.4f}, see-radius {rv:.4f} D={rv*G:.3f} Gpc")
fig,ax=plt.subplots(figsize=(13,11));th=np.linspace(0,2*np.pi,721)
P=np.r_[phi[ok],-phi[ok][::-1]];RE=np.r_[re[ok],re[ok][::-1]];RT=np.r_[rt[ok],rt[ok][::-1]]
ax.plot(RE*np.cos(P),RE*np.sin(P),'b-',lw=2.5,label=f'where the CMB we see was emitted (at emission time); from us: 0° {de[0]*G:.2f}, 90° {de[90]*G:.2f}, 180° {de[180]*G:.2f} Gpc')
ax.plot(RT*np.cos(P),RT*np.sin(P),'r-',lw=2.5,label=f'same matter today, moved straight out from CoE; from us: 0° {dt[0]*G:.2f}, 90° {dt[90]*G:.2f}, 180° {dt[180]*G:.2f} Gpc')
ax.plot(xc+rv*np.cos(th),rv*np.sin(th),'g--',lw=2,label=f'how far we can see: circle of radius c×(time since we crossed) = {rv*G:.2f} Gpc')
ax.plot(R*np.cos(th),R*np.sin(th),'k:',lw=1.5,label=f'recombination radius {R} D ({R*G:.2f} Gpc)')
for s in (1,-1):
    for i in range(0,181,15):
        if ok[i]: ax.plot([re[i]*np.cos(s*phi[i]),rt[i]*np.cos(s*phi[i])],[re[i]*np.sin(s*phi[i]),rt[i]*np.sin(s*phi[i])],color='gray',lw=.8)
ax.plot(xo,0,'k*',ms=16);ax.annotate('us today',(xo,0),(8,8),textcoords='offset points',fontsize=11)
ax.plot(xc,0,'ko',ms=8,mfc='none');ax.annotate('us when we crossed',(xc,0),(8,-16),textcoords='offset points',fontsize=10)
ax.plot(0,0,'k+',ms=18,mew=2);ax.annotate('CoE',(0,0),(8,8),textcoords='offset points',fontsize=11)
X=np.r_[RE*np.cos(P),RT*np.cos(P),xc-rv,xc+rv,xo,0];Y=np.r_[RE*np.sin(P),RT*np.sin(P),rv,-rv];m=0.08*(X.max()-X.min());ax.set_xlim(X.min()-m,X.max()+m);ax.set_ylim(Y.min()-m,Y.max()+m);ax.set_aspect('equal');ax.grid(alpha=.3);ax.set_xlabel('distance along CoE–us line (D = 4.58 Gpc)');ax.set_ylabel('distance (D)')
xl=ax.get_xlim();ax.text(xl[1],0,'0° anti-CoE side →',ha='right',va='bottom',fontsize=12);ax.text(xl[0],-0.02,'← 180° CoE side',ha='left',va='top',fontsize=12)
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
fig.savefig(f'surface_{tag}.png',dpi=110,bbox_inches='tight')
