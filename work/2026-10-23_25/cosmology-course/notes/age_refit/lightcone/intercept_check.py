"""Check figure for Chris's photon-intercept construction (letter eqs 19-23; same as notes/cmb-cancel.py).
Origin O at (0,0); observer 1 at (1,0) at emission; emitters on the concentric circle of radius q around observer 1.
Frozen radial speeds v=k r (beta1 = observer speed, c=1), optional gravity mu via v^2 = beta1^2 r0^2 + 2 mu (1/r - 1/r0) (s=1 here).
Photon: straight line source->observer at reception, (v_o^2-1) t^2 + 2 dx v_o t + |d|^2 = 0, positive root t_f."""
import sys,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
beta1,q,mu=[float(v) for v in sys.argv[1:4]] if len(sys.argv)>3 else (0.6,0.4,0.0)
th=np.radians(np.arange(0,360,30.)); xe=1+q*np.cos(th); ye=q*np.sin(th); re=np.hypot(xe,ye)
vo2=beta1**2; ve2=beta1**2*re**2+2*mu*(1/re-1/re); vo,ve=np.sqrt(vo2),np.sqrt(ve2)
dx,dy=1-xe,-ye; a=vo2-1; b=2*dx*vo; c=dx**2+dy**2; tf=(-b-np.sqrt(b*b-4*a*c))/(2*a)
nx,ny=(dx+vo*tf)/tf,dy/tf; bex,bey=ve*xe/re,ve*ye/re
Dm=(1/np.sqrt(1-vo2))*(1-vo*nx)/((1/np.sqrt(1-ve2))*(1-(bex*nx+bey*ny))); Dg=np.sqrt((1-2*mu/re)/(1-2*mu)); S=Dm*Dg
fig,ax=plt.subplots(figsize=(19,12)); ax.set_aspect('equal'); cc=np.linspace(0,2*np.pi,300)
ax.plot(1+q*np.cos(cc),q*np.sin(cc),c='orange',lw=4,label=f'emitters: concentric circle radius q={q} around observer 1 (at emission)')
T=tf.max(); ax.plot([0,1+vo*T],[0,0],'k-',lw=1); ax.plot([1,1+vo*T],[0,0],'b-',lw=5,label=f'observer 1 worldline, v_o={vo:.2f}c')
ax.plot(1,0,'o',c='b',ms=16)
for i in range(len(th)):
    ax.plot([xe[i],xe[i]+ve[i]*xe[i]/re[i]*tf[i]],[ye[i],ye[i]+ve[i]*ye[i]/re[i]*tf[i]],'-',c='green',lw=2)   # emitter worldline during flight
    ax.annotate('',xy=(1+vo*tf[i],0),xytext=(xe[i],ye[i]),arrowprops=dict(arrowstyle='->',color='C3',lw=2.5))
    ax.plot(xe[i],ye[i],'o',c='orange',ms=12); ax.plot(1+vo*tf[i],0,'x',c='C3',ms=14,mew=3)
    ax.annotate(f'{np.degrees(th[i]):.0f}°: t_f={tf[i]:.3f}\nS={S[i]:.3f}',(xe[i],ye[i]),fontsize=13,xytext=(8,8 if ye[i]>=0 else -36),textcoords='offset points')
ax.plot(0,0,'kx',ms=20,mew=4); ax.text(0.02,-0.08,'O (origin / CoE)',fontsize=18)
ax.set_title(f'Exact photon-intercept paths (c=1): green = emitter worldline during flight, red = photon path c·t_f to observer at reception (×)\nβ1={beta1}, q={q}, μ={mu}; θ=0 points away from O. S = D_motion·D_grav (eqs 19–23)',fontsize=18)
ax.set_xlabel('x / d1',fontsize=20); ax.set_ylabel('y / d1',fontsize=20); ax.tick_params(labelsize=18); ax.grid(); ax.legend(fontsize=14,loc='upper left')
fig.tight_layout(); fig.savefig('intercept_check.png',dpi=80)
print('t_f',np.round(tf,4).tolist()); print('S',np.round(S,4).tolist())
