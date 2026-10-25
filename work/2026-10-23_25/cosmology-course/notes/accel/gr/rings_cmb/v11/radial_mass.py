import sys,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
u=5.0;f=float(sys.argv[1]);tag=sys.argv[2];Rrec=0.5;t0=1/u;tc=Rrec/u;G=4.2/(t0-tc);Gy=G*3.26156
M=f*0.5*u**2   # same mass fraction of critical: M/(u^2 D/2) at our distance D, D=1,c=1 (v11 base: 0.00853 = 0.0154*0.5*1.052^2)
a=np.radians(np.arange(0,360,2.0));rho=0.4*Rrec;xc,yc=u*tc+rho*np.cos(a),rho*np.sin(a)
r0=np.hypot(xc,yc);ph=np.arctan2(yc,xc);v0=r0/tc
r0=np.r_[r0,u*tc];v0=np.r_[v0,u]   # last = us
ts=[tc,tc+(t0-tc)/3,tc+2*(t0-tc)/3,t0];n=len(r0)
s=solve_ivp(lambda t,y:np.r_[y[n:],-M/y[:n]**2],(tc,t0),np.r_[r0,v0],t_eval=ts,rtol=1e-11,atol=1e-13)
R=s.y[:n];V=s.y[n:]
print(f'{tag}: fraction {f}, M={M:.5f} (D=1,c=1), D={G:.2f} Gpc, today {t0*Gy:.3f} Gyr; our speed crossing {u} c -> today {V[-1,-1]:.4f} c; min particle speed today {V[:-1,-1].min():.4f} c')
th=np.linspace(0,2*np.pi,721);fig,ax=plt.subplots(figsize=(14,10))
cols=['lightcoral','indianred','firebrick','darkred']
for j,t in enumerate(ts):
    x,y=R[:-1,j]*np.cos(ph),R[:-1,j]*np.sin(ph);us=R[-1,j]
    ax.plot(np.r_[x,x[0]]*G,np.r_[y,y[0]]*G,'-',color=cols[j],lw=2);ax.plot(us*G,0,'k*',ms=11)
    ax.annotate(f'{t*Gy:.1f} Gyr',(x[45]*G,y[45]*G),(0,6),textcoords='offset points',color=cols[j],fontsize=9,ha='center')
    rc_=us-R[90,j];ra=R[0,j]-us;rp=np.hypot(x[45]-us,y[45])
    print(f' t={t*Gy:6.2f} Gyr: from us: CoE side {rc_*G:.3f}, anti-CoE side {ra*G:.3f}, 90deg {rp*G:.3f} Gpc; CoE/anti -1 = {100*(rc_/ra-1):+.3f}%; length/width -1 = {100*((rc_+ra)/(2*np.hypot(x,y-0).max()*0+2*np.abs(y).max())-1):+.3f}%')
for i in range(0,180,10): ax.plot([0,R[i,-1]*np.cos(ph[i])*G],[0,R[i,-1]*np.sin(ph[i])*G],color='gray',lw=.5)
ax.plot(Rrec*G*np.cos(ph),Rrec*G*np.sin(ph),'o',color='mediumblue',ms=3,zorder=5)
ax.fill(Rrec*G*np.cos(th),Rrec*G*np.sin(th),color='lightblue',alpha=.3);ax.plot(Rrec*G*np.cos(th),Rrec*G*np.sin(th),'b-',lw=1.2)
ust=R[-1,-1];gr=t0-tc;ax.plot((ust+gr*np.cos(th))*G,gr*np.sin(th)*G,'g-',lw=2)
dmin=np.hypot(R[:-1,-1]*np.cos(ph)-ust,R[:-1,-1]*np.sin(ph)).min();print(f' green {gr*G:.3f} Gpc vs closest red today {dmin*G:.3f} Gpc -> inside {dmin>=gr}; closest red to CoE today {R[:-1,-1].min()*G:.3f} vs Rrec {Rrec*G:.3f} Gpc')
ax.plot(0,0,'k+',ms=18,mew=2);ax.annotate('CoE / central mass',(0,0),(-40,10),textcoords='offset points')
ax.plot([],[],'-',color='darkred',lw=2,label=f'red ring of particles at {", ".join(f"{t*Gy:.1f}" for t in ts)} Gyr (light→dark; star = us), each moving only along its gray radial line')
ax.plot([],[],'g-',lw=2,label=f'how far we can see today: c × time since we crossed = {gr*G:.2f} Gpc, centred on us')
ax.plot([],[],'o',color='mediumblue',label='dark blue: where each particle path crossed the recombination sphere')
ax.plot([],[],'b-',label=f'recombination sphere 0.5 D = {Rrec*G:.1f} Gpc (shaded = still plasma)')
ax.set_aspect('equal');ax.grid(alpha=.3);ax.set_xlabel('Gpc (anti-CoE side to the right)');ax.set_ylabel('Gpc')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
ax.set_title(f'Central mass {f:g} of critical; our speed at crossing 5 c; D = {G:.0f} Gpc',fontsize=11,pad=-14)
fig.savefig(f'red_radial_mass_{tag}.png',dpi=110,bbox_inches='tight')
