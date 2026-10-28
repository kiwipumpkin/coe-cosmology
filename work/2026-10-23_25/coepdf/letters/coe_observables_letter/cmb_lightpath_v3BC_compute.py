# v3B and v3C: Newtonian CoE model (same fit/lock/per-shell recombination as cmb_lightpath_v2_compute.py), two light rules.
# Shared: shell trajectories r(q,t), v(q,t) from the prof_core acceleration/RK4-log-t scheme (coasting x=q t, v=q before taui);
#   advected (time-dependent) mass M(<r,t)=M_lag(q), M_lag(q)=MT(r(q,t0)); Phi(r,t)=-mu[M(<r,t)/r + int_r^inf dM/r'];
#   t_rec(q) where rho(q,t)=rho_us,now*(3000/2.725)^3  (T∝1/a, rho∝a^-3).  Units: lengths D, times t0, c=0.024359 D/t0.
# B: photon velocity = v_flow(x,t) + c k, k=-n (n = look direction), fixed; integrate backward from us (shell q=1 today) in ln t
#    until t = t_rec(q(x,t)).  ln(1+z_dop) = int (k.grad v.k) dt  (exp of integral of local Doppler between neighbouring shells),
#    observer comoving with the flow.  grav = (Phi_e-Phi_o)/c^2 ; ISW = (2/c^2) int dPhi/dt|_x dt along the ray (GR factor 2).
#    T_obs = T_rec * exp(-ln(1+z_dop)) * (1 + grav + ISW); T_rec fixed (universal density threshold) -> no T_rec-variation term.
# C: emitter = shell whose PRESENT distance from us along n is c(t0 - t_rec(q)); event (r(q,t_rec) along that shell's radial line,
#    t_rec); linear Doppler -n.(V_e(t_rec)-V_o(t0))/c (flagged); grav as B; ISW along the straight segment emission->us with time
#    linear in path fraction (not a null path; flagged); T_rec fixed -> no extra term.
import os,sys,json,numpy as np,importlib.util as iu
from scipy.integrate import cumulative_trapezoid
BUB='/workspace/cosmology/papers/spreadsheet_model/bubble';sys.path.insert(0,BUB);cwd=os.getcwd();os.chdir(BUB)
import prof_core as pc
from bubble_core import RKM,GPC,TH0
sp=iu.spec_from_file_location('cw','/workspace/cosmology-course/notes/catwise-eta-density-migkas.py');cw=iu.module_from_spec(sp);sp.loader.exec_module(cw);os.chdir(cwd)
FIT=json.load(open('/workspace/coepdf/joint_solve/phase2/p1_migkas_D34.96.json'))['34.96'];Dg,eta,mc=FIT['D'],FIT['eta'],FIT['mc'];RV=float(RKM/GPC)
L=cw.L_from_delta_eta(Dg,RV,eta,cw.DELTA_USE);rr=pc.run('eta',dict(DR=Dg/RV,LD=L/Dg,eta=eta),mc,0.0)
mu,MT,S,Ht=rr['mu'],rr['MT'],rr['S'],rr['Ht'];t0=Ht*TH0;cG=299792.458*3.15576e16/3.0857e22;c=1/(Dg/t0/cG)
T0,TR=2.725,3000.;TAUI,NST=0.01,400
q=np.geomspace(1e-4,4.0,4000)
def a_r(r):r=np.maximum(r,1e-12);return -mu*np.interp(r,pc.RG,MT)/r**2
s0=np.log(TAUI);h=-s0/NST;x=q*TAUI;v=q.copy();s=s0;XR=[x];VV=[v];tt=[TAUI]
for n in range(NST):
    t1,t2,t4=np.exp(s),np.exp(s+h/2),np.exp(s+h)
    k1=t1*v;l1=t1*a_r(x);k2=t2*(v+h/2*l1);l2=t2*a_r(x+h/2*k1);k3=t2*(v+h/2*l2);l3=t2*a_r(x+h/2*k2);k4=t4*(v+h*l3);l4=t4*a_r(x+h*k3)
    x=x+h/6*(k1+2*k2+2*k3+k4);v=v+h/6*(l1+2*l2+2*l3+l4);s+=h;XR.append(x);VV.append(v);tt.append(np.exp(s))
XR=np.array(XR);VV=np.array(VV);lt=np.log(np.array(tt))
r0=XR[-1];Mlag=np.interp(r0,pc.RG,MT);dM=np.gradient(Mlag,q)
Mq=4*np.pi*np.interp(r0,pc.RG,S)*r0**2*np.gradient(r0,q)
def state(t):
    """r(q),v(q) at time t (arrays over q)"""
    if t<=TAUI:return q*t,q.copy()
    u=np.clip((np.log(t)-lt[0])/h,0,NST-1e-9);i=int(u);w=u-i
    return (1-w)*XR[i]+w*XR[i+1],(1-w)*VV[i]+w*VV[i+1]
def phiq(rq):
    tail=cumulative_trapezoid((dM/rq)[::-1],q[::-1],initial=0)[::-1]*-1  # int_q^inf dM/r
    return -mu*(Mlag/rq+tail)
def field(t,r):
    rq,vq=state(t);qq=np.interp(r,rq,q);vv=np.interp(r,rq,vq);dv=np.interp(r,rq,np.gradient(vq,rq));ph=np.interp(r,rq,phiq(rq))
    return qq,vv,dv,ph
def dphidt(t,r,e=1e-3):return (field(t*(1+e),r)[3]-field(t*(1-e),r)[3])/(2*e*t)
rho_now_us=np.interp(np.interp(1.0,q,r0),pc.RG,S);rho_rec=rho_now_us*(TR/T0)**3
def rhoq(t):rq,_=state(t);return Mq/(4*np.pi*rq**2*np.gradient(rq,q))
# t_rec(q): coasting regime check
R0=rhoq(TAUI);trec=TAUI*(R0/rho_rec)**(1/3)   # coasting: rho∝t^-3
late=np.where(R0>=rho_rec)[0];RR=np.array([Mq/(4*np.pi*XR[i]**2*np.gradient(XR[i],q)) for i in range(NST+1)])
for jq in late:
    kk=np.where(RR[:,jq]<rho_rec)[0]
    if len(kk)==0:trec[jq]=np.nan;continue
    kk=kk[0];w=(np.log(RR[kk-1,jq])-np.log(rho_rec))/(np.log(RR[kk-1,jq])-np.log(RR[kk,jq]));trec[jq]=np.exp(lt[kk-1]+w*h)
print('late shells',len(late),'never',int(np.isnan(trec).sum()),'trec range',np.nanmin(trec),np.nanmax(trec))
okq=np.isfinite(trec);trec=np.interp(q,q[okq],trec[okq])  # shells that never reach threshold (if any) get edge value; count printed
ro=float(np.interp(1.0,q,r0));vo=float(np.interp(1.0,q,VV[-1]));Xo=np.array([ro,0.]);Phio=field(1.0,np.array([ro]))[3][0]
psi=np.radians(np.linspace(0,180,181))
# ---------- B ----------
NB=1500;sB=np.linspace(0,np.log(1e-4),NB+1);resB=[]
for p in psi:
    n=np.array([np.cos(p),np.sin(p)]);k=-n;X=Xo.copy();lnz=0.;isw=0.;hit=None
    for j in range(NB):
        ta,tb=np.exp(sB[j]),np.exp(sB[j+1]);dt=tb-ta     # negative
        def vel(X,t):
            r=np.hypot(*X);qq,vv,dv,ph=field(t,np.array([r]));rh=X/r;return vv[0]*rh+c*k,(dv[0]*(k@rh)**2+vv[0]/r*(1-(k@rh)**2)),dphidt(t,np.array([r]))[0],qq[0]
        u1,g1,i1,q1=vel(X,ta);Xm=X+0.5*dt*u1;u2,g2,i2,_=vel(Xm,0.5*(ta+tb));u3,g3,i3,_=vel(X+0.5*dt*u2,0.5*(ta+tb));u4,g4,i4,q4=vel(X+dt*u3,tb)
        Xn=X+dt/6*(u1+2*u2+2*u3+u4);lnz+=-dt/6*(g1+2*g2+2*g3+g4);isw+=-dt/6*(i1+2*i2+2*i3+i4)
        f1=ta-np.interp(q1,q,trec);f2=tb-np.interp(q4,q,trec)
        if f2<=0:
            w=f1/(f1-f2);te=ta+w*dt;Xe=X+w*(Xn-X);lnz-=(1-w)*(-dt/6*(g1+2*g2+2*g3+g4));isw-=(1-w)*(-dt/6*(i1+2*i2+2*i3+i4))
            re=np.hypot(*Xe);qe,ve,_,phe=field(te,np.array([re]));hit=(te,Xe,qe[0],phe[0]);break
        X=Xn
    if hit is None:resB.append([np.degrees(p)]+[np.nan]*8);continue
    te,Xe,qe,phe=hit;grav=(phe-Phio)/c**2;iswT=2*isw/c**2
    resB.append([np.degrees(p),te,Xe[0],Xe[1],qe,lnz,grav,iswT,np.hypot(*(Xe-Xo))])
resB=np.array(resB)
# ---------- C ----------
resC=[]
for p in psi:
    n=np.array([np.cos(p),np.sin(p)]);sv=c*(1-np.interp(1.0,q,trec))
    for it in range(50):
        P=Xo+n*sv;qe=np.interp(np.hypot(*P),r0,q);te=np.interp(qe,q,trec);sn=c*(1-te)
        if abs(sn-sv)<1e-14:break
        sv=sn
    rh=P/np.hypot(*P);rq,vq=state(te);re=np.interp(qe,q,rq);ve=np.interp(qe,q,vq);Xe=rh*re;Ve=rh*ve
    Vo=np.array([vo,0.]);kin=-n@(Ve-Vo)/c
    phe=field(te,np.array([re]))[3][0];grav=(phe-Phio)/c**2
    lam=np.linspace(0,1,801);iswI=0.
    ts=te+(1-te)*lam;pts=Xe[None]+(Xo-Xe)[None]*lam[:,None];dpt=np.array([dphidt(ts[i],np.array([np.hypot(*pts[i])]))[0] for i in range(len(lam))])
    iswT=2*np.trapezoid(dpt,ts)/c**2
    resC.append([np.degrees(p),te,Xe[0],Xe[1],qe,kin,grav,iswT,np.hypot(*(Xe-Xo)),sv])
resC=np.array(resC)
np.savez('cmb_lightpath_v3BC.npz',resB=resB,resC=resC,cols_B='psi,t_e,xe,ye,q_e,ln1pz_dop,grav,ISW,dist_emission_event_from_us_now',
  cols_C='psi,t_e,xe,ye,q_e,kin_lin,grav,ISW,dist_emission_event,s_label',c=c,D=Dg,t0=t0,ro=ro,vo=vo,q=q,trec=trec)
print('B nan:',int(np.isnan(resB[:,1]).sum()),'lnz range',np.nanmin(resB[:,5]),np.nanmax(resB[:,5]),'grav',np.nanmin(resB[:,6]),np.nanmax(resB[:,6]),'isw',np.nanmin(resB[:,7]),np.nanmax(resB[:,7]))
print('C kin',resC[:,5].min(),resC[:,5].max(),'grav',resC[:,6].min(),resC[:,6].max(),'isw',resC[:,7].min(),resC[:,7].max())
