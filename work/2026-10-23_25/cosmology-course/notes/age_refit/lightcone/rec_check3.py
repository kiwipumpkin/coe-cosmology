"""rec_check3: Chris's fixed-R_rec construction with light at c RELATIVE TO LOCAL MATTER (c=1, lengths in R_rec, CoE fixed at origin).
Flow: every element obeys dr/ds=U, dU/ds=-mu/r^2 after crossing R_rec (U=U_rec there) -> steady radial flow u(r)=sqrt(U_rec^2-2mu(1/R_rec-1/r)).
Photon in moving medium (Galilean ray optics, Hamiltonian H=c|k|+k.u): dx/dt=c k/|k|+u, dk/dt=-grad(k.u). Speed c relative to local matter.
Steady flow -> H conserved; local frequency nu'=c|k| -> expansion redshift 1+z_exp=|k_e|/|k_o| (from flow stretching only).
Motion part: peculiar velocities set to 0 (all matter on the flow) -> 1+z_mot=1. Gravity: weak-field exp(mu/r_e... ) see below.
Emission event: backward ray from us first reaches r=R_rec (element crossing there emits once, 3000 K)."""
import numpy as np, json
from scipy.integrate import solve_ivp
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
Rrec,Urec,mu,ts=1.0,3.0,1.0,2.0
u=lambda r: np.sqrt(Urec**2-2*mu*(1/Rrec-1/r)); du=lambda r: -mu/(r*r*u(r))
our=solve_ivp(lambda s,y:[y[1],-mu/y[0]**2],(0,ts),[Rrec,Urec],rtol=1e-12,atol=1e-12); D=our.y[0,-1]; Uo=our.y[1,-1]
def rhs(t,y):
    x=y[:2];k=y[2:];r=np.hypot(*x);kx=k@x/r;U=u(r)
    return np.r_[k/np.linalg.norm(k)+U*x/r, -(du(r)*kx*x/r+U*(k/r-kx*x/r**2))]
hit=lambda t,y: np.hypot(*y[:2])-Rrec; hit.terminal=True
rows=[];paths=[]
for th in np.arange(0,360,30.):
    n=np.array([np.cos(np.radians(th)),np.sin(np.radians(th))])
    s=solve_ivp(rhs,(ts,ts-60),np.r_[D,0,-n],events=hit,rtol=1e-10,atol=1e-12,max_step=0.01)
    if not s.t_events[0].size: rows.append(dict(theta=th,hit=False));continue
    te=s.t_events[0][0];E=s.y_events[0][0];ke=np.linalg.norm(E[2:])
    zexp=ke-1; zgr=np.exp(mu/Rrec-mu/D)-1          # weak-field potential climb Rrec->D (Schwarzschild sqrt invalid: 2mu/Rrec=2)
    ex=s.y[0];ey=s.y[1];L=np.sum(np.hypot(np.diff(ex),np.diff(ey)))
    rows.append(dict(theta=th,hit=True,t_e=te,dt=ts-te,local_path=ts-te,CoE_frame_path=L,emit_xy=E[:2].tolist(),emit_angle_CoE=float(np.degrees(np.arctan2(E[1],E[0]))%360),
        beyond_CoE=bool(E[0]<0),z_exp=zexp,z_mot=0.0,z_grav=zgr,z_tot=(1+zexp)*(1+zgr)-1,T_obs=3000/((1+zexp)*(1+zgr))))
    paths.append((th,ex,ey))
fig,ax=plt.subplots(figsize=(20,13));ax.set_aspect('equal');c_=np.linspace(0,2*np.pi,400)
ax.add_patch(plt.Circle((0,0),Rrec,color='#ffe0e0'));ax.plot(np.cos(c_),np.sin(c_),'r--',lw=3,label='R_rec (fixed); inside = plasma')
R=ts;ax.plot(D+R*np.cos(c_),R*np.sin(c_),c='orange',lw=3,label=f'visible circle R=c t*={R} (local)')
for rr,l in ((D-R,'D−R'),(D+R,'D+R')): ax.plot(rr*np.cos(c_),rr*np.sin(c_),':',c='gray',lw=2);ax.text(rr*.71,rr*.71,l,fontsize=18)
ax.plot([Rrec,D],[0,0],'b-',lw=4,label=f'our path (U {Urec}c → {Uo:.2f}c)');ax.plot(D,0,'b*',ms=28,label=f'us now D={D:.2f}')
for (th,ex,ey),q in zip(paths,[q for q in rows if q['hit']]):
    ax.plot(ex,ey,'g-',lw=2.5);ax.plot(ex[-1],ey[-1],'go',ms=12)
    ax.annotate(f"{th:.0f}°: t_e={q['t_e']:.2f}, z={q['z_tot']:.2f}",(ex[-1],ey[-1]),fontsize=14,xytext=(8,8 if ey[-1]>=0 else -18),textcoords='offset points')
ax.plot(0,0,'k*',ms=24);ax.text(.05,-.3,'CoE',fontsize=20)
nb=sum(q.get('beyond_CoE',False) for q in rows)
ax.set_title(f'rec_check3: emission once at R_rec crossing, light at c relative to local flow (green = photon path in CoE frame)\nU_rec=3c, μ=1, t*=2. Beyond CoE: {nb}/12',fontsize=19)
ax.tick_params(labelsize=16);ax.grid();ax.legend(fontsize=15,loc='lower left');fig.tight_layout();fig.savefig('rec_check3.png',dpi=80)
json.dump(dict(D=D,U_now=Uo,rows=rows),open('rec_check3.json','w'),indent=1,default=float)
print(D,Uo)
for q in rows: print({k:(round(v,3) if isinstance(v,float) else v) for k,v in q.items() if k!='emit_xy'})
