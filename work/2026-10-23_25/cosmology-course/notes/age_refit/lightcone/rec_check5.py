"""rec_check5: rec_check4 + scale factor a(t)=D(t)/D(0) (D(0)=R_rec=1) for the co-moving region. Light at c in local proper distance:
comoving radius chi(t)=int_0^t dt'/a(t'); proper radius a(t)chi(t). Emitters at comoving offset chi(t*) n, emitted at t=0 (region crossing).
1+z=(a_obs/a_emit)(1+z_g)(1+z_m); z_g, z_m as rec_check4 at emitter's present proper position D+a chi n."""
import numpy as np,json
from scipy.integrate import solve_ivp,quad
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
Rrec,Urec,mu,ts=1.,3.,1.,2.
s=solve_ivp(lambda t,y:[y[1],-mu/y[0]**2],(0,ts),[Rrec,Urec],dense_output=True,rtol=1e-12,atol=1e-12)
D=lambda t:s.sol(t)[0];Uf=lambda t:s.sol(t)[1];a=lambda t:D(t)/Rrec
chi=lambda t:quad(lambda x:1/a(x),0,t)[0]
Do,Uo,ao=D(ts),Uf(ts),a(ts);X=chi(ts);Rp=ao*X
u=lambda r:np.sqrt(Urec**2-2*mu*(1/Rrec-1/r))
rows=[]
for th in np.arange(0,360,30.):
    n=np.array([np.cos(np.radians(th)),np.sin(np.radians(th))]);E=np.array([Do,0])+Rp*n;re=np.linalg.norm(E)
    zg=np.sqrt((1-2*mu/Do)/(1-2*mu/re))-1;vp=(u(re)-Uo)*E/re;zm=(1+vp@n)/np.sqrt(1-vp@vp)-1;za=ao-1
    zt=(1+za)*(1+zg)*(1+zm)-1
    rows.append(dict(theta=th,r_emit_now=re,z_exp=za,z_grav=zg,z_mot=zm,z_tot=zt,T_obs=3000/(1+zt)))
T=np.array([q['T_obs'] for q in rows]);th=np.radians(np.arange(0,360,30.))
d=lambda arr:(arr-arr.mean())/arr.mean();dT=d(T)
c1=2*np.mean(dT*np.cos(th));c2=2*np.mean(dT*np.cos(2*th))
info=dict(D=Do,U=Uo,a_obs=ao,chi=X,R_proper=Rp,ct=ts,Tmean=T.mean(),stretch_needed=3000/2.725,extra_factor=T.mean()/2.725,dipole=c1,quad=c2,rows=rows)
c_=np.linspace(0,2*np.pi,400);steps=[0.5,1.0,1.5,2.0];cols=['C0','C1','C2','C3']
fig,axs=plt.subplots(4,1,figsize=(18,26))
for ax,t,cl in zip(axs,steps,cols):
    dd=D(t);R=a(t)*chi(t);ax.set_aspect('equal');ax.add_patch(plt.Circle((0,0),Rrec,color='#ffe0e0'));ax.plot(np.cos(c_),np.sin(c_),'r--',lw=2)
    ax.add_patch(Wedge((0,0),dd+R,0,360,width=2*R,color=cl,alpha=.12))
    for rr in (dd-R,dd+R): ax.plot(rr*np.cos(c_),rr*np.sin(c_),':',c=cl,lw=2.5)
    ax.plot(dd+R*np.cos(c_),R*np.sin(c_),c=cl,lw=3)
    for t_ in np.arange(0,360,30.):
        p=np.array([dd,0])+R*np.array([np.cos(np.radians(t_)),np.sin(np.radians(t_))]);ax.annotate('',xy=(dd,0),xytext=p,arrowprops=dict(arrowstyle='->',color='k',lw=1.6));ax.plot(*p,'o',c=cl,ms=9)
    ax.plot(dd,0,'b*',ms=24);ax.plot(0,0,'k*',ms=20);ax.text(.1,-.5,'CoE',fontsize=20)
    ax.set_title(f't*={t}: D={dd:.2f}, a={a(t):.2f}, U={Uf(t):.2f}c; proper R={R:.2f} (c t*={t}), comoving χ={chi(t):.3f}',fontsize=20)
    ax.set_xlim(-11,12);ax.set_ylim(-4.5,4.5);ax.tick_params(labelsize=16);ax.grid()
fig.suptitle('rec_check5: stretching co-moving region a(t)=D(t)/R_rec; photons at c in proper distance (straight in comoving coords)',fontsize=22)
fig.tight_layout(rect=[0,0,1,.97]);fig.savefig('rec_check5.png',dpi=70)
json.dump(info,open('rec_check5.json','w'),indent=1,default=float)
for k,v in info.items():
    if k!='rows':print(k,round(float(v),4))
for q in rows:print({k:round(float(v),4) for k,v in q.items()})
# --- steering: explicit travel time and clock-rate split
for q in rows:
    th_=np.radians(q['theta']);n=np.array([np.cos(th_),np.sin(th_)]);E=np.array([Do,0])+Rp*n;re=np.linalg.norm(E)
    vp=(u(re)-Uo)*E/re;g=1/np.sqrt(1-vp@vp)
    q.update(t_emit=0.0,t_travel=ts,clock_grav=1+q['z_grav'],clock_kin_gamma=g,doppler_los=1+vp@n,stretch=ao)
json.dump(info,open('rec_check5.json','w'),indent=1,default=float)
for q in rows:print(q['theta'],round(q['clock_grav'],4),round(q['clock_kin_gamma'],6),round(q['doppler_los'],4))
