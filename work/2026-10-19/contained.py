import sys,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
npz,out,R,lab=sys.argv[1],sys.argv[2],float(sys.argv[3]),sys.argv[4]
d=np.load(npz);o=d['o'];re,rt,phi=d['re'],d['rt'],d['phi'];xo,xc,t0,tc=[float(d[k]) for k in('xo','xc','t0','tc')]
rv=t0-tc;G=4.2/rv;Gy=G*3.26156;ok=np.isfinite(o[:,1])
print(f'{lab}: found {ok.sum()}/{len(ok)}; circle {rv:.5f} D; D={G:.4f} Gpc; today {t0*Gy:.3f} Gyr; crossing {tc*Gy:.3f} Gyr')
dd=lambda r,p:np.hypot(r*np.cos(p)-xo,r*np.sin(p));de,dt=dd(re,phi),dd(rt,phi)
for i in(0,90,180): print(f'{i}: emit dist {de[i]*G:.2f}, today dist {dt[i]*G:.2f}, light travel {o[i,3]*G:.2f} Gpc, emit t {o[i,2]*Gy:.2f} Gyr, 1+z {o[i,1]:.3f}')
P=np.r_[phi[ok],-phi[ok][::-1]];RE=np.r_[re[ok],re[ok][::-1]]*G;RT=np.r_[rt[ok],rt[ok][::-1]]*G;X0,XC=xo*G,xc*G
DD=np.hypot(RT*np.cos(P)-X0,RT*np.sin(P));print(f'min red dist from us {DD.min():.4f} Gpc vs circle 4.2000 -> contained: {DD.min()>=4.2}')
x,y=RT*np.cos(P),RT*np.sin(P);A=np.c_[2*x,2*y,np.ones_like(x)];c=np.linalg.lstsq(A,x*x+y*y,rcond=None)[0];rc=np.sqrt(c[2]+c[0]**2+c[1]**2);dev=np.abs(np.hypot(x-c[0],y-c[1])-rc)
print(f'red circle fit: centre ({c[0]:.3f},{c[1]:.3f}) Gpc, radius {rc:.3f} Gpc, max deviation {dev.max():.3f} Gpc = {100*dev.max()/rc:.2f}%')
fig,ax=plt.subplots(figsize=(13,11));th=np.linspace(0,2*np.pi,721)
ax.plot(RE*np.cos(P),RE*np.sin(P),'b-',lw=2.5,label=f'where the CMB we see was emitted (at emission time); from us: 0° {de[0]*G:.2f}, 90° {de[90]*G:.2f}, 180° {de[180]*G:.2f} Gpc')
ax.plot(x,y,'r-',lw=2.5,label=f'same matter today, moved straight out from CoE; from us: 0° {dt[0]*G:.2f}, 90° {dt[90]*G:.2f}, 180° {dt[180]*G:.2f} Gpc')
ax.plot(X0+4.2*np.cos(th),4.2*np.sin(th),'g--',lw=2,label='How far we can see: c × time since we crossed = 4.20 Gpc, centred on us today')
ax.plot(R*G*np.cos(th),R*G*np.sin(th),'k:',lw=1.5);ax.annotate(f'recombination radius {R} D = {R*G:.2f} Gpc (dotted)',(0,-R*G),(-40,-16),textcoords='offset points',fontsize=10)
for s in(1,-1):
    for i in range(0,181,15):
        if not ok[i]: continue
        xr,yr=rt[i]*G*np.cos(s*phi[i]),rt[i]*G*np.sin(s*phi[i])
        ax.plot([re[i]*G*np.cos(s*phi[i]),xr],[re[i]*G*np.sin(s*phi[i]),yr],color='gray',lw=.8);ax.plot([xr,X0],[yr,0],color='darkorange',lw=1)
        if s>0 and i in(0,90,180):
            L=np.hypot(xr-X0,yr);print(f'ray {i}: length {L:.2f} Gpc, crosses green at '+(f'({X0+4.2*(xr-X0)/L:.2f},{4.2*yr/L:.2f})' if L>4.2 else 'none'))
ax.plot([],[],color='darkorange',lw=1,label='straight light ray from red point to us today (every 15°; gray = radial lines from CoE)')
ax.plot(X0,0,'k*',ms=16);ax.annotate('us today',(X0,0),(8,8),textcoords='offset points',fontsize=11)
ax.plot(XC,0,'ko',ms=8,mfc='none');ax.annotate('us when we crossed',(XC,0),(8,-16),textcoords='offset points',fontsize=10)
ax.plot(0,0,'k+',ms=18,mew=2);ax.annotate('CoE',(0,0),(8,8),textcoords='offset points',fontsize=11)
X=np.r_[RE*np.cos(P),x,X0-4.2,X0+4.2,-R*G];Y=np.r_[RE*np.sin(P),y,4.2,-4.2,-R*G];m=0.08*(X.max()-X.min())
ax.set_xlim(X.min()-m,X.max()+m);ax.set_ylim(Y.min()-m,Y.max()+m);ax.set_aspect('equal');ax.grid(alpha=.3)
ax.set_xlabel('Gpc along the CoE–us line');ax.set_ylabel('Gpc')
xl=ax.get_xlim();ax.text(xl[1],0,'0° anti-CoE side →',ha='right',va='bottom',fontsize=12);ax.text(xl[0],-0.15,'← 180° CoE side',ha='left',va='top',fontsize=12)
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
ax.set_title(f'{lab}; D = {G:.2f} Gpc (set so the model seeing radius is 4.2 Gpc); recombination radius {R} D',fontsize=11,pad=-14)
fig.savefig(out,dpi=110,bbox_inches='tight')
