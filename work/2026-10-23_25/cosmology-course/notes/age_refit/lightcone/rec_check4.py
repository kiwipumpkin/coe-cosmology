"""rec_check4: co-moving region between CoE-centred shells D-R, D+R; light straight at c in that frame. c=1, R_rec=1, U_rec=3, mu=1, t*=2.
Our D(t): dr/dt=U, dU/dt=-mu/r^2 from (R_rec,U_rec) at t=0 (we crossed at t=0). Emitters recombined at t=0 at frame offset R n (R=c t*), at rest in frame.
Redshift: (a) gravity D_g: 1+z_g=sqrt((1-2mu/r_o)/(1-2mu/r_e)) [Schwarzschild static clocks], r_e=|(D+R cos th, R sin th)| (emitter's CoE radius, evaluated now);
(b) peculiar: flow speed u(r)=sqrt(U_rec^2-2mu(1/R_rec-1/r)); v_pec=(u(r_e)-U_o) radial at emitter (only the radial-speed DIFFERENCE; frame's own angular spread excluded);
1+z_m=gamma(1+v_pec.n), n=unit vector us->emitter (relativistic Doppler, emitter receding along n -> redshift)."""
import numpy as np,json
from scipy.integrate import solve_ivp
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
Rrec,Urec,mu,ts=1.,3.,1.,2.
s=solve_ivp(lambda t,y:[y[1],-mu/y[0]**2],(0,ts),[Rrec,Urec],dense_output=True,rtol=1e-12,atol=1e-12)
D=lambda t:s.sol(t)[0];Uf=lambda t:s.sol(t)[1];Do,Uo=D(ts),Uf(ts)
u=lambda r:np.sqrt(Urec**2-2*mu*(1/Rrec-1/r))
rows=[]
for th in np.arange(0,360,30.):
    n=np.array([np.cos(np.radians(th)),np.sin(np.radians(th))]);E=np.array([Do,0])+ts*n;re=np.linalg.norm(E)
    zg=np.sqrt((1-2*mu/Do)/(1-2*mu/re))-1
    vp=(u(re)-Uo)*E/re;vl=vp@n;zm=(1+vl)/np.sqrt(1-vp@vp)-1
    rows.append(dict(theta=th,t_emit=0.0,path=ts,r_emit=re,v_pec=float(np.linalg.norm(vp)*np.sign(u(re)-Uo)),v_los=vl,z_grav=zg,z_mot=zm,z_tot=(1+zg)*(1+zm)-1,T_obs=3000/((1+zg)*(1+zm)),beyond_CoE=bool(E[0]<0)))
c_=np.linspace(0,2*np.pi,400);steps=[0.5,1.0,1.5,2.0];cols=['C0','C1','C2','C3']
fig,axs=plt.subplots(len(steps),1,figsize=(18,26))
for ax,t,cl in zip(axs,steps,cols):
    d,R=D(t),t;ax.set_aspect('equal');ax.add_patch(plt.Circle((0,0),Rrec,color='#ffe0e0'));ax.plot(np.cos(c_),np.sin(c_),'r--',lw=2)
    from matplotlib.patches import Wedge
    ax.add_patch(Wedge((0,0),d+R,0,360,width=2*R,color=cl,alpha=.12))
    for rr in (d-R,d+R): ax.plot(rr*np.cos(c_),rr*np.sin(c_),':',c=cl,lw=2.5)
    ax.plot(d+R*np.cos(c_),R*np.sin(c_),c=cl,lw=3)
    for th in np.arange(0,360,30.):
        p=np.array([d,0])+R*np.array([np.cos(np.radians(th)),np.sin(np.radians(th))]);ax.annotate('',xy=(d,0),xytext=p,arrowprops=dict(arrowstyle='->',color='k',lw=1.6));ax.plot(*p,'o',c=cl,ms=9)
    ax.plot(d,0,'b*',ms=24);ax.plot(0,0,'k*',ms=20);ax.text(.1,-.5,'CoE',fontsize=20)
    ax.set_title(f't*={t}: D={d:.2f}, U={Uf(t):.2f}c, visible R=c t*={R}; region = shells D−R={d-R:.2f} to D+R={d+R:.2f}',fontsize=20)
    ax.set_xlim(-9,9.5);ax.set_ylim(-3,3);ax.tick_params(labelsize=16);ax.grid()
fig.suptitle('rec_check4: co-moving region (shaded) carried outward; straight photons at c inside it from emitters on our circle to us',fontsize=22)
fig.tight_layout(rect=[0,0,1,.97]);fig.savefig('rec_check4.png',dpi=70)
json.dump(dict(D=Do,U=Uo,rows=rows),open('rec_check4.json','w'),indent=1,default=float)
print(Do,Uo,D(0.5),D(1),D(1.5));[print({k:round(float(v),4) for k,v in q.items()}) for q in rows]
