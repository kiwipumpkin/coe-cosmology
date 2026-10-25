import numpy as np,json
src=open('steps.py').read(); exec(src[:src.index('import matplotlib')])
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
Dt=lambda t: r_of(1.0,t)
def rhs(t,y,n): r=np.hypot(*y); return (v_of_r(r,t)/r)*y + C*n
ev=lambda t,y,n: np.hypot(*y)-RREC; ev.terminal=True
targets=[5,10,15,20,25,Dt(1.0)*DG]; out=[]; rays={}
for DGt in targets:
    tA=1.0 if DGt==targets[-1] else brentq(lambda t: Dt(t)*DG-DGt,trec*1.01,1.0)
    L=[]
    for beta in np.arange(0,360,15):
        n=-np.array([np.cos(np.radians(beta)),np.sin(np.radians(beta))])
        sb=solve_ivp(rhs,(tA,T[5]),[Dt(tA),0.0],args=(n,),events=ev,rtol=1e-10,atol=1e-13,max_step=2e-3)
        te=sb.t_events[0][0]; xe=sb.y_events[0][0]
        sf=solve_ivp(rhs,(te,tA),xe,args=(n,),rtol=1e-11,atol=1e-14,dense_output=True,max_step=2e-3)
        err=np.hypot(*(sf.y[:,-1]-[Dt(tA),0]))*DG; assert err<1e-4,(DGt,beta,err)
        tt=np.linspace(te,tA,600); Y=sf.sol(tt); Yr=Y-np.vstack([[Dt(x) for x in tt],0*tt])
        L.append(dict(te=te,xe=xe,Y=Y,Yr=Yr,err=err))
    ang=np.degrees([np.arctan2(r['xe'][1],r['xe'][0]) for r in L]); tes=np.array([r['te'] for r in L])
    span=ang.max()-ang.min(); rays[DGt]=(tA,L)
    out.append(dict(D_Gpc=float(Dt(tA)*DG),t_Gyr=tA*t0,tstar_Gyr=(tA-trec)*t0,R_ct_Gpc=C*(tA-trec)*DG,
        te_min_Gyr=tes.min()*t0,te_max_Gyr=tes.max()*t0,ang_span_deg=float(span),patch_Gpc=float(np.radians(span)*RREC*DG),max_err_Gpc=max(r['err'] for r in L)))
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig=plt.figure(figsize=(26,38)); gs=fig.add_gridspec(4,3,height_ratios=[1.0,1.1,1,1])
A=fig.add_subplot(gs[0,:]); cols=plt.cm.tab10(np.arange(6))
ph=np.linspace(0,2*np.pi,600); A.plot(RREC*DG*np.cos(ph),RREC*DG*np.sin(ph),color='purple',lw=4,label='recombination sphere R_rec = 2 Gpc')
A.plot(0,0,'*',color='r',ms=30,label='CoE')
for (DGt,(tA,L)),c,o in zip(rays.items(),cols,out):
    for r in L: A.plot(r['Y'][0]*DG,r['Y'][1]*DG,'-',color=c,lw=1.5,alpha=.8)
    E=np.array([r['xe'] for r in L])*DG; A.plot(E[:,0],E[:,1],'o',color=c,ms=9,mec='k')
    A.plot(Dt(tA)*DG,0,'*',color=c,ms=26,mec='k',label=f"D = {o['D_Gpc']:.2f} Gpc (patch {o['ang_span_deg']:.1f}°)")
Z=fig.add_subplot(gs[1,:])
for (DGt,(tA,L)),c,o in zip(rays.items(),cols,out):
    E=np.array([r['xe'] for r in L]); ang=np.degrees(np.arctan2(E[:,1],E[:,0])); te=np.array([r['te'] for r in L])*t0
    Z.plot(ang,te,'o-',color=c,ms=11,mec='k',lw=2,label=f"D = {o['D_Gpc']:.1f}")
Z.axhline(trec*t0,color='k',ls=':',lw=2,label='we cross R_rec')
Z.set_xlabel('emission point angle on R_rec, seen from CoE (deg; 0 = our direction)',fontsize=24); Z.set_ylabel('emission time t_e (Gyr)',fontsize=24)
Z.tick_params(labelsize=20); Z.legend(fontsize=18,ncol=4); Z.set_title('Which part of the recombination sphere feeds us at each D',fontsize=28)
A.set_aspect('equal'); A.set_xlim(-3,29); A.set_ylim(-3,3)
A.set_xlabel('x from CoE (Gpc)',fontsize=26); A.set_ylabel('y (Gpc)',fontsize=26); A.tick_params(labelsize=22)
A.set_title('CoE frame, true scale: rays arriving at us at 6 different D',fontsize=30); A.legend(fontsize=19,loc='upper right',ncol=2)
for k,((DGt,(tA,L)),c,o) in enumerate(zip(rays.items(),cols,out)):
    B=fig.add_subplot(gs[2+k//3,k%3]); R=o['R_ct_Gpc']
    for r in L: B.plot(r['Yr'][0]*DG,r['Yr'][1]*DG,'-',color=c,lw=2)
    B.plot(*(np.array([r['Yr'][:,0] for r in L]).T*DG),'o',color=c,ms=7,mec='k')
    B.plot(R*np.cos(ph),R*np.sin(ph),'--',color='green',lw=3); B.plot(0,0,'*',color='k',ms=22)
    B.set_aspect('equal'); B.set_xlim(-1.15*R,1.15*R); B.set_ylim(-1.15*R,1.15*R); B.tick_params(labelsize=17)
    B.set_title(f"D={o['D_Gpc']:.1f}, t={o['t_Gyr']:.2f} Gyr\nR=ct*={R:.3f} Gpc (green)",fontsize=22)
    B.set_xlabel('Gpc rel. to us',fontsize=19)
fig.tight_layout(); fig.savefig('rays_vs_D.png',dpi=65)
json.dump(out,open('rays_vs_D.json','w'),indent=1)
for o in out: print({k:round(v,4) for k,v in o.items()})
