import sys,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt; from ic2 import *
beta1,q,mu,U=(0.3,0.4,0.0,1.5)
ell=np.arange(0,360,30.); r=S_ic2(beta1,q,mu,1.0,U,ell)
fig,ax=plt.subplots(figsize=(19,12)); ax.set_aspect('equal'); cc=np.linspace(0,2*np.pi,300)
ax.plot(1+q*np.cos(cc),q*np.sin(cc),c='orange',lw=4,label=f'emitters at emission: circle q={q} around us (physical = comoving at t=0)')
T=r['tf'].max(); tt=np.linspace(0,T,200); tau=lambda t: np.log1p(U*t)/U; a=lambda t:1+U*t
ax.plot(a(tt)*(1+r['vo']*tau(tt)),0*tt,'b-',lw=5,label=f'our worldline (global U={U}c + peculiar {beta1}c)')
for i in range(len(ell)):
    t=np.linspace(0,r['tf'][i],100); ta=tau(t); xe,ye=r['xe0'][i],r['ye0'][i]; re0=np.hypot(xe,ye)
    ax.plot(a(t)*(xe+r['ve'][i]*xe/re0*ta),a(t)*(ye+r['ve'][i]*ye/re0*ta),'g-',lw=2)               # emitter worldline during flight
    xo=1+r['vo']*r['tau'][i]; f=ta/r['tau'][i]; X=a(t)*(xe+(xo-xe)*f); Y=a(t)*(ye*(1-f))               # photon: straight in comoving, stretched by a(t)
    ax.plot(X,Y,'-',c='C3',lw=2.5); ax.plot(X[-1],Y[-1],'x',c='C3',ms=14,mew=3); ax.plot(xe,ye,'o',c='orange',ms=12)
    ax.annotate(f"{ell[i]:.0f}°: t_f={r['tf'][i]:.2f}, path {r['tf'][i]:.2f} c\n1+z_U={1+r['zU'][i]:.2f}, S={r['S'][i]:.3f}",(xe,ye),fontsize=12,xytext=(8,8 if ye>=0 else -36),textcoords='offset points')
ax.plot(0,0,'kx',ms=20,mew=4); ax.text(0.02,-0.1,'O (CoE)',fontsize=18)
ax.set_title(f'Intercept v2: shell recedes globally at U={U}c (carries us + emitters); peculiar speeds <c enter D_motion\nred = photon path as flown (straight in shell coords, stretched by expansion) to us at reception (×); green = emitter worldlines',fontsize=18)
ax.set_xlabel('x / d1',fontsize=20); ax.set_ylabel('y / d1',fontsize=20); ax.tick_params(labelsize=18); ax.grid(); ax.legend(fontsize=14,loc='upper left')
fig.tight_layout(); fig.savefig('intercept_check2.png',dpi=80)
print('tf',np.round(r['tf'],3).tolist()); print('S',np.round(r['S'],4).tolist()); print('Dm',np.round(r['Dm'],4).tolist())
