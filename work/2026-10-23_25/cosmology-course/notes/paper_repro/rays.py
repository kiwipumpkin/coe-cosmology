import numpy as np,json
src=open('steps.py').read(); exec(src[:src.index('import matplotlib')])
from scipy.integrate import solve_ivp
tN=1.0; Dt=lambda t: r_of(1.0,t); Ut=lambda t: v_of_r(Dt(t),t)
def rhs(t,y,n):
    r=np.hypot(*y); return (v_of_r(r,t)/r)*y + C*n
def back(beta):          # arrival direction beta (deg, photon comes FROM this sky direction); n = -(cos,sin)
    n=-np.array([np.cos(np.radians(beta)),np.sin(np.radians(beta))])
    ev=lambda t,y,n: np.hypot(*y)-RREC; ev.terminal=True
    s=solve_ivp(rhs,(tN,T[5]),[Dt(tN),0.0],args=(n,),events=ev,rtol=1e-10,atol=1e-13,dense_output=True,max_step=2e-3)
    return n,s.t_events[0][0],s.y_events[0][0],s
rays=[];rows=[]
for beta in np.arange(0,360,30):
    n,te,xe,sb=back(beta)
    sf=solve_ivp(rhs,(te,tN),xe,args=(n,),rtol=1e-11,atol=1e-14,dense_output=True,max_step=2e-3)   # forward shot
    err=np.hypot(*(sf.y[:,-1]-[Dt(tN),0]))*DG; assert err<1e-4,(beta,err)
    tt=np.linspace(te,tN,1500); Y=sf.sol(tt); Dtt=np.array([Dt(x) for x in tt])
    Yr=Y-np.vstack([Dtt,0*Dtt])
    L_coe=np.sum(np.hypot(*np.diff(Y)))*DG; L_rest=np.sum(np.hypot(*np.diff(Yr)))*DG
    qe=q_of_r(RREC,te)
    rays.append(dict(beta=beta,tt=tt,Y=Y,Yr=Yr,xe=xe,qe=qe,te=te))
    rows.append(dict(beta=int(beta),t_e_Gyr=te*t0,travel_Gyr=(tN-te)*t0,c_travel_Gpc=C*(tN-te)*DG,L_CoE_Gpc=L_coe,L_rest_Gpc=L_rest,arrival_err_Gpc=err,
        emit_x=xe[0]*DG,emit_y=xe[1]*DG))
Rstar=C*(tN-trec)*DG
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
fig,(A,B)=plt.subplots(2,1,figsize=(20,24))
norm=plt.Normalize(trec*t0*0.95,t0)
def cl(ax,xy,tt):
    pts=xy.T.reshape(-1,1,2); seg=np.concatenate([pts[:-1],pts[1:]],1)
    lc=LineCollection(seg,cmap='plasma',norm=norm,lw=3); lc.set_array(tt[:-1]*t0); ax.add_collection(lc); return lc
tg=np.linspace(trec*0.9,tN,400); DD=np.array([Dt(x) for x in tg])
A.plot(DD*DG,0*DD,'-',color='orange',lw=5,label='our track D(t)')
for r in rays:
    tt=np.linspace(r['te'],tN,300); rr_=np.array([r_of(r['qe'],x) for x in tt]); u=r['xe']/np.hypot(*r['xe'])
    A.plot(rr_*u[0]*DG,rr_*u[1]*DG,'-',color='gray',lw=1.5,alpha=.7)
    lc=cl(A,r['Y']*DG,r['tt']); A.plot(*(r['xe']*DG),'o',color='cyan',ms=10,mec='k')
A.plot([],[],'-',color='gray',label="emitters' later tracks (radial)"); A.plot([],[],'o',color='cyan',mec='k',label='emission point (on R_rec)')
A.plot(Dt(tN)*DG,0,'*',color='lime',ms=30,mec='k',label=f'arrival: us at D(now) = {Dt(tN)*DG:.2f} Gpc')
th=np.linspace(-0.2,0.2,100); A.plot(RREC*DG*np.cos(th),RREC*DG*np.sin(th),'-',color='purple',lw=3,label='R_rec = 2 Gpc')
A.set_xlim(0,29); A.set_ylim(-1.0,1.0); A.set_xlabel('x from CoE (Gpc)',fontsize=24); A.set_ylabel('y (Gpc)  [stretched ×~14]',fontsize=24)
A.set_title('CoE frame: traced photons (velocity = u(r,t) + c·n̂)',fontsize=28); A.legend(fontsize=19,loc='upper left'); A.tick_params(labelsize=20)
cb=fig.colorbar(lc,ax=A); cb.set_label('time along ray (Gyr)',fontsize=22); cb.ax.tick_params(labelsize=18)
for r in rays: lc=cl(B,r['Yr']*DG,r['tt']); B.plot(*(r['Yr'][:,0]*DG),'o',color='cyan',ms=10,mec='k')
B.plot(0,0,'*',color='lime',ms=34,mec='k'); ph=np.linspace(0,2*np.pi,400)
B.plot(Rstar*np.cos(ph),Rstar*np.sin(ph),'--',color='green',lw=3,label=f'c·t* circle, R = {Rstar:.3f} Gpc')
B.set_aspect('equal'); m=Rstar*1.1; B.set_xlim(-m,m); B.set_ylim(-m,m)
B.set_xlabel('x relative to us (Gpc)',fontsize=24); B.set_ylabel('y relative to us (Gpc)',fontsize=24); B.tick_params(labelsize=20)
B.set_title('Our rest frame: same rays converging on us (equal axes)',fontsize=28); B.legend(fontsize=20,loc='upper right')
cb=fig.colorbar(lc,ax=B); cb.set_label('time along ray (Gyr)',fontsize=22); cb.ax.tick_params(labelsize=18)
fig.tight_layout(); fig.savefig('traced_rays.png',dpi=70)
for x in rows: print({k:(round(v,5) if isinstance(v,float) else v) for k,v in x.items()})
print('ct* circle',Rstar,'trec',trec*t0)
json.dump(dict(rows=rows,Rstar_Gpc=Rstar),open('traced_rays.json','w'),indent=1)
