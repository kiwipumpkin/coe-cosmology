import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.load('surface_mass_R05_42circle.npz');o=d['o'];re,rt,phi=d['re'],d['rt'],d['phi'];xo,xc,t0,tc=[float(d[k]) for k in('xo','xc','t0','tc')]
rv=t0-tc;G=4.2/rv;Gy=G*3.26156;R=0.5
print(f'D={G:.4f} Gpc, seeing radius {rv:.5f} D, D/c={Gy:.4f} Gyr, today t={t0*Gy:.3f} Gyr, crossing t={tc*Gy:.3f} Gyr, since crossing {rv*Gy:.3f} Gyr; old D 4.58 -> today {t0*4.58*3.26156:.3f} Gyr')
ok=np.isfinite(o[:,1]);dd=lambda r,p:np.hypot(r*np.cos(p)-xo,r*np.sin(p));de,dt=dd(re,phi),dd(rt,phi)
for i in(0,90,180): print(f'{i}: emit dist {de[i]*G:.2f} Gpc, today dist {dt[i]*G:.2f} Gpc, light travel {o[i,3]*G:.2f} Gpc, emit t {o[i,2]*Gy:.2f} Gyr ({o[i,2]:.4f} D/c), 1+z {o[i,1]:.3f}')
P=np.r_[phi[ok],-phi[ok][::-1]];RE=np.r_[re[ok],re[ok][::-1]]*G;RT=np.r_[rt[ok],rt[ok][::-1]]*G;X0,XC=xo*G,xc*G
fig,ax=plt.subplots(figsize=(13,11));th=np.linspace(0,2*np.pi,721)
ax.plot(RE*np.cos(P),RE*np.sin(P),'b-',lw=2.5,label=f'where the CMB we see was emitted (at emission time); from us: 0° {de[0]*G:.2f}, 90° {de[90]*G:.2f}, 180° {de[180]*G:.2f} Gpc')
ax.plot(RT*np.cos(P),RT*np.sin(P),'r-',lw=2.5,label=f'same matter today, moved straight out from CoE; from us: 0° {dt[0]*G:.2f}, 90° {dt[90]*G:.2f}, 180° {dt[180]*G:.2f} Gpc')
ax.plot(X0+4.2*np.cos(th),4.2*np.sin(th),'g--',lw=2,label='How far we can see, 4.2 Gpc')
ax.plot(R*G*np.cos(th),R*G*np.sin(th),'k:',lw=1.5,label=f'recombination radius 0.5 D = {R*G:.2f} Gpc')
for s in(1,-1):
    for i in range(0,181,15): ax.plot([RE[i]*0+re[i]*G*np.cos(s*phi[i]),rt[i]*G*np.cos(s*phi[i])],[re[i]*G*np.sin(s*phi[i]),rt[i]*G*np.sin(s*phi[i])],color='gray',lw=.8)
ax.plot(X0,0,'k*',ms=16);ax.annotate('us today',(X0,0),(8,8),textcoords='offset points',fontsize=11)
ax.plot(XC,0,'ko',ms=8,mfc='none');ax.annotate('us when we crossed',(XC,0),(8,-16),textcoords='offset points',fontsize=10)
ax.plot(0,0,'k+',ms=18,mew=2);ax.annotate('CoE',(0,0),(8,8),textcoords='offset points',fontsize=11)
X=np.r_[RE*np.cos(P),RT*np.cos(P),X0-4.2,X0+4.2,0];Y=np.r_[RE*np.sin(P),RT*np.sin(P),rv*G,-rv*G];m=0.08*(X.max()-X.min())
ax.set_xlim(X.min()-m,X.max()+m);ax.set_ylim(Y.min()-m,Y.max()+m);ax.set_aspect('equal');ax.grid(alpha=.3)
ax.set_xlabel('Gpc along the CoE–us line');ax.set_ylabel('Gpc')
xl=ax.get_xlim();ax.text(xl[1],0,'0° anti-CoE side →',ha='right',va='bottom',fontsize=12);ax.text(xl[0],-0.15,'← 180° CoE side',ha='left',va='top',fontsize=12)
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
ax.set_title(f'CoE-to-us distance D = {G:.2f} Gpc (rescaled so the seeing radius is 4.2 Gpc); with core mass, recombination radius 0.5 D',fontsize=11,y=1.0,pad=-14)
fig.savefig('surface_mass_R05_Dscaled_v2.png',dpi=110,bbox_inches='tight')
