import json,os,numpy as np
N='/workspace/cosmology-course/notes/'
src=open(N+'catwise-density-levelset-migkas.py').read(); ns={'__name__':'x'}; exec(src[:src.index('def main')],ns)
os.chdir(N+'paper_repro'); pc=ns['pc']; import bubble_core as bc
p=json.load(open(N+'catwise-density-levelset-fit.json'))['migkas_fits']['soft_fixed15.2']
a=ns['_pc_args']('soft',p['DR'],p['LD'],p['shape']); rr=pc.run('soft',a,p['mc'],0.0)
mu,MT,RG,S=rr['mu'],rr['MT'],pc.RG,rr['S']
DG=p['D_Gpc']; t0=rr['Ht']*bc.TH0          # Gyr: unit time of the code
C=0.306601*t0/DG                            # c in code units (D/t0); c=0.306601 Gpc/Gyr
# radial version of prof_core.integrate (same ICs x=q*taui, v=q at t=taui, same log-time RK4, same acc with mub=0)
taui=0.01; NS=20000; s0=np.log(taui); h=-s0/NS; T=np.exp(s0+h*np.arange(NS+1))
def acc(r): return -mu*np.interp(r,RG,MT)/r**2
def traj(q):
    q=np.atleast_1d(np.asarray(q,float)); x=q*taui; v=q.copy(); X=[x.copy()]; Vv=[v.copy()]; s=s0
    for n in range(NS):
        t1,t2,t4=np.exp(s),np.exp(s+.5*h),np.exp(s+h)
        l1=t1*acc(x);k1=t1*v; l2=t2*acc(x+.5*h*k1);k2=t2*(v+.5*h*l1)
        l3=t2*acc(x+.5*h*k2);k3=t2*(v+.5*h*l2); l4=t4*acc(x+h*k3);k4=t4*(v+h*l3)
        x=x+h/6*(k1+2*k2+2*k3+k4); v=v+h/6*(l1+2*l2+2*l3+l4); s+=h; X.append(x.copy()); Vv.append(v.copy())
    return np.array(X),np.array(Vv)
# family of layers q in [qmin,qmax]; positions r(q,t) monotone in q
QG=np.linspace(0.02,1.6,1581); XR,VR=traj(QG)
Xc,Vc,_=pc.integrate(np.array([[1.0,0.0]]),mu,RG,MT,0.0,rr['Rb'],taui,200)   # code's own observer
i1=np.argmin(abs(QG-1)); assert abs(XR[-1,i1]-Xc[0,0])<1e-5 and abs(VR[-1,i1]-Vc[0,0])<1e-5,(XR[-1,i1],Xc[0,0])
def at(t,arr): j=np.searchsorted(T,t); w=(t-T[j-1])/(T[j]-T[j-1]); return arr[j-1]*(1-w)+arr[j]*w
def r_of(q,t): return np.interp(q,QG,at(t,XR))
def v_of_r(r,t): return np.interp(r,at(t,XR),at(t,VR))
def q_of_r(r,t): return np.interp(r,at(t,XR),QG)
RREC=2.0/DG                                   # chosen: R_rec = 2 Gpc
Dus=XR[:,i1]; trec=np.interp(RREC,Dus,T)      # time we crossed R_rec
tk=trec+(1-trec)*np.arange(1,5)/4; dt=tk[1]-tk[0]
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig,axs=plt.subplots(4,2,figsize=(30,46),gridspec_kw=dict(width_ratios=[1,1]))
Z=1.15*C*(1-trec)*DG                           # zoom half-width (Gpc), same every panel
SM=6.0/ (VR[0,i1]/C)  ; SZ=0.35*Z               # arrow scales: main Gpc per c ; zoom Gpc per c
g=np.linspace(-1.25*DG,1.25*DG,500); xx,yy=np.meshgrid(g,g); dens=np.interp(np.hypot(xx,yy)/DG,RG,S)/np.interp(1,RG,S)
prev=None; log=[]; th=np.linspace(0,2*np.pi,721)
for k,(t,(A,B)) in enumerate(zip(tk,axs)):
    D=r_of(1.0,t); U=v_of_r(D,t); R=C*(t-trec); qi,qo=q_of_r(D-R,t),q_of_r(D+R,t)
    vals=dict(t_Gyr=t*t0,D_Gpc=D*DG,U_c=U/C,R_Gpc=R*DG,dR_Gpc=(C*dt*DG if k else R*DG))
    for ax in (A,B):
        ax.pcolormesh(g,g,dens,cmap='Greys',alpha=.35,shading='auto',vmin=dens.min(),vmax=min(dens.max(),3))
        ax.plot(0,0,'*',color='r',ms=30); ax.plot(RREC*DG*np.cos(th),RREC*DG*np.sin(th),'-',color='purple',lw=3)
        for rr_,c_ in ((D-R,'tab:blue'),(D+R,'tab:red')): ax.plot(rr_*DG*np.cos(th),rr_*DG*np.sin(th),'-',color=c_,lw=3)
        if prev is not None:
            for qq,c_ in zip(prev,('tab:blue','tab:red')):
                rp=r_of(qq,t); ax.plot(rp*DG*np.cos(th),rp*DG*np.sin(th),'--',color=c_,lw=3)
        ax.plot(D*DG+R*DG*np.cos(th),R*DG*np.sin(th),'-',color='green',lw=4)
        ax.plot(D*DG,0,'o',color='orange',ms=18,mec='k',zorder=5)
    # main panel: global velocities
    for rr_,c_,lab in ((D,'orange','U (us)'),(D-R,'tab:blue','inner shell'),(D+R,'tab:red','outer shell')):
        y0={'U (us)':0,'inner shell':-0.12*DG,'outer shell':0.12*DG}[lab]
        A.annotate('',xy=(rr_*DG+SM*v_of_r(rr_,t)/C,y0),xytext=(rr_*DG,y0),arrowprops=dict(arrowstyle='-|>',lw=4,color=c_))
    A.set_xlim(-0.3*DG,1.25*DG);A.set_ylim(-0.6*DG,0.6*DG);A.set_aspect('equal')
    A.text(-0.28*DG,0.5*DG,f"t = {vals['t_Gyr']:.3f} Gyr   D = {vals['D_Gpc']:.2f} Gpc   U/c = {vals['U_c']:.1f}\nR = c t* = {vals['R_Gpc']:.4f} Gpc   ΔR = {vals['dR_Gpc']:.4f} Gpc",fontsize=24,va='top',bbox=dict(fc='w',ec='k'))
    A.set_xlabel('x (Gpc)',fontsize=24);A.set_ylabel('y (Gpc)',fontsize=24);A.tick_params(labelsize=20)
    A.text(0.3*DG,-0.55*DG,f'arrows: {SM:.3f} Gpc per c',fontsize=20)
    # zoom: peculiar velocities on visible circle
    ang=np.radians(np.arange(0,360,30)); P=np.c_[D+R*np.cos(ang),R*np.sin(ang)]; rp=np.hypot(*P.T)
    Vp=(v_of_r(rp,t)/rp)[:,None]*P - np.array([U,0.0])
    if R>0:
        for (x,y),(vx,vy) in zip(P,Vp): B.annotate('',xy=(x*DG+SZ*vx/C,y*DG+SZ*vy/C),xytext=(x*DG,y*DG),arrowprops=dict(arrowstyle='-|>',lw=3,color='k'))
    B.set_xlim(D*DG-Z,D*DG+Z);B.set_ylim(-Z,Z);B.set_aspect('equal')
    B.set_xlabel('x (Gpc)',fontsize=24);B.tick_params(labelsize=20)
    B.set_title(f'zoom around us (±{Z:.2f} Gpc); black: emitter velocity relative to us, {SZ:.3f} Gpc per c',fontsize=20)
    log.append(dict(vals,qi=qi,qo=qo,Vp_max_c=float(np.max(np.hypot(*Vp.T))/C)))
    # assertions: plotted = computed
    assert np.isclose(D*DG,vals['D_Gpc']) and np.isclose(np.hypot(*(P[0]-[D,0])),R) and np.isclose(r_of(qi,t),D-R,atol=1e-6) and np.isclose(r_of(qo,t),D+R,atol=1e-6)
    if k: assert np.isclose(vals['R_Gpc']-log[k-1]['R_Gpc'],C*dt*DG)
    prev=(qi,qo)
# CMB emitters on last panel: each point on today's circle, traced back along its radial path to R_rec
t=tk[-1]; D=r_of(1.0,t); R=C*(t-trec); ang=np.radians(np.arange(0,360,30)); P=np.c_[D+R*np.cos(ang),R*np.sin(ang)]
rp=np.hypot(*P.T); qe=q_of_r(rp,t); te=np.array([np.interp(RREC,XR[:,np.argmin(abs(QG-q))],T) for q in qe]); E=RREC*P/rp[:,None]
for ax in axs[-1]:
    ax.plot(P[:,0]*DG,P[:,1]*DG,'D',color='magenta',ms=12,mec='k'); ax.plot(E[:,0]*DG,E[:,1]*DG,'D',color='magenta',ms=12,mec='k')
    for a_,e_ in zip(P,E): ax.plot([a_[0]*DG,e_[0]*DG],[a_[1]*DG,e_[1]*DG],':',color='magenta',lw=1.5)
axs[-1,0].text(-0.28*DG,-0.35*DG,f"CMB emitters (magenta): crossed R_rec at t = {te.min()*t0:.4f}–{te.max()*t0:.4f} Gyr\n(we crossed at {trec*t0:.4f} Gyr)",fontsize=20,bbox=dict(fc='w'))

Dte=np.array([r_of(1.0,x) for x in te]); A0=axs[-1,0]
for e_,dte,tt in zip(E,Dte,te):
    A0.plot([e_[0]*DG,D*DG],[e_[1]*DG,0],'-',color='gold',lw=2,alpha=.9)
    A0.plot(dte*DG,0,'o',color='orange',ms=10,alpha=.35,mec='k')
A0.plot([],[],'-',color='gold',lw=2,label='light path: emitter at emission → us now')
A0.plot([],[],'o',color='orange',alpha=.35,label='us at each emission time')
A0.text(0.3*DG,0.3*DG,f"travel time t−t_e = {(t-te.max())*t0:.3f}–{(t-te.min())*t0:.3f} Gyr\nc·(travel) = {C*(t-te.max())*DG:.3f}–{C*(t-te.min())*DG:.3f} Gpc",fontsize=20,bbox=dict(fc='w'))
A0.legend(fontsize=18,loc='lower right')
assert np.allclose(np.hypot(*E.T),RREC)
from matplotlib.lines import Line2D
H=[Line2D([],[],marker='*',color='r',ls='',ms=20,label='CoE'),Line2D([],[],color='purple',lw=3,label=f'R_rec = {RREC*DG:.1f} Gpc'),
   Line2D([],[],marker='o',color='orange',ls='',ms=15,label='us at D(t)'),Line2D([],[],color='green',lw=4,label='visible circle R = c t*'),
   Line2D([],[],color='tab:blue',lw=3,label='shell D−R (now)'),Line2D([],[],color='tab:red',lw=3,label='shell D+R (now)'),
   Line2D([],[],color='gray',ls='--',lw=3,label='previous step\'s layers, integrated'),Line2D([],[],marker='D',color='magenta',ls='',ms=12,label='CMB emitters / crossing points')]
fig.legend(handles=H,loc='upper center',ncol=4,fontsize=22,bbox_to_anchor=(0.5,0.975))
fig.suptitle('Soft-core best fit, prof_core dynamics: 4 steps (same axes every row)',fontsize=30,y=0.995)
fig.tight_layout(rect=(0,0,1,0.95)); fig.savefig(N+'paper_repro/steps_fitted.png',dpi=70)
print(json.dumps(dict(t0_Gyr=t0,c_code=C,trec_Gyr=trec*t0,R_rec_Gpc=RREC*DG,D_today_int=r_of(1,1)*DG,steps=log,te_Gyr=list(te*t0),D_at_te_Gpc=list(Dte*DG)),indent=1))
