import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
u=10**0.021927885054513485   # our speed in c (v11 best fit), no mass: r = u*t
Rrec=0.5;t0=1/u;tc=Rrec/u;k=1.1   # red radius at crossing = k * c * t_cross (chosen)
G=4.2/(t0-tc);Gy=G*3.26156
print(f'u={u:.4f}c t0={t0:.4f} tc={tc:.4f} D/c; D={G:.4f} Gpc; age today {t0*Gy:.3f} Gyr; crossing {tc*Gy:.3f} Gyr')
th=np.linspace(0,2*np.pi,721);eps=k*tc
ts=[tc,tc+(t0-tc)/3,tc+2*(t0-tc)/3,t0];sty=[':','-.','--','-'];al=[.4,.6,.8,1]
fig,ax=plt.subplots(figsize=(13,11))
for j,(t,s,a) in enumerate(zip(ts,sty,al)):
    xs,ys=u*tc+eps*np.cos(th),eps*np.sin(th)           # red matter at crossing
    xr,yr=xs*t/tc,ys*t/tc                                # each point moved radially from CoE, r proportional to t
    us=u*t;fr=np.hypot(xr-us,yr);rg=t-tc
    cx,cy=xr[:-1].mean(),yr[:-1].mean();dev=np.abs(np.hypot(xr-cx,yr-cy)-fr.mean()).max()
    rate=eps/tc  # d(red radius)/dt in c
    print(f't={t*Gy:6.3f} Gyr: red radius {fr.mean()*G:.3f} Gpc (centre offset from us {np.hypot(cx-us,cy)*G:.1e}), max dev from circle {dev*G:.1e} Gpc; green {rg*G:.3f} Gpc; green inside: {fr.min()>=rg}; red growth {rate:.3f} c, green 1 c')
    ax.plot(xr*G,yr*G,'r',ls=s,lw=2.2,alpha=a);ax.plot((us+rg*np.cos(th))*G,rg*np.sin(th)*G,'g',ls=s,lw=2,alpha=a)
    ax.plot(us*G,0,'k*',ms=12,alpha=a);ax.annotate(f'us, {t*Gy:.1f} Gyr',(us*G,0),(-20,-18-14*j),textcoords='offset points',fontsize=9)
ax.plot(0,0,'k+',ms=18,mew=2);ax.annotate('CoE',(0,0),(8,8),textcoords='offset points',fontsize=11)
ax.plot(Rrec*G*np.cos(th),Rrec*G*np.sin(th),color='gray',lw=1)
ax.plot([],[],'r-',lw=2.2,label=f'red: circle of matter around us at crossing (radius 1.1 × c × crossing time), carried out by expansion; grows at {eps/tc:.2f} c')
ax.plot([],[],'g-',lw=2,label='green: how far we can see, c × time since we crossed, centred on us; grows at 1 c')
ax.plot([],[],'gray',lw=1,label=f'recombination radius 0.5 D = {Rrec*G:.2f} Gpc')
ax.plot([],[],'k',ls=':',label=f'line style = time: dotted {ts[0]*Gy:.1f}, dash-dot {ts[1]*Gy:.1f}, dashed {ts[2]*Gy:.1f}, solid today {ts[3]*Gy:.1f} Gyr')
ax.set_aspect('equal');ax.grid(alpha=.3);ax.set_xlabel('Gpc (CoE at 0; anti-CoE side to the right)');ax.set_ylabel('Gpc')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
ax.set_title(f'No mass; D = {G:.2f} Gpc, today = {t0*Gy:.1f} Gyr',fontsize=11,pad=-14)
fig.savefig('nomass_red_green_simple.png',dpi=110,bbox_inches='tight')
