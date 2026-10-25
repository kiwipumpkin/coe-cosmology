import numpy as np,json,sys
src=open('steps.py').read(); src=src[:src.index('import matplotlib')]
src=src.replace("taui=0.01; NS=20000","taui=1e-6; NS=80000").replace("i1=np.argmin(abs(QG-1)); assert","i1=np.argmin(abs(QG-1)); DIFF=(XR[-1,i1]-Xc[0,0],VR[-1,i1]-Vc[0,0]); 0 and")
exec(src)
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
Dt=lambda t: r_of(1.0,t); dq=2e-3
def rhs(t,y,n): r=np.hypot(*y); return (v_of_r(r,t)/r)*y + C*n
def lay(q,t): r=r_of(q,t); return r,(r_of(q+dq,t)-r_of(q-dq,t))
def phi(r): return float(pc.phi_coe(np.atleast_1d(r),mu,MT)[0])
def rays(Rr,step):
    ev=lambda t,y,n: np.hypot(*y)-Rr; ev.terminal=True; out=[]
    for beta in np.arange(0,360,step):
        n=-np.array([np.cos(np.radians(beta)),np.sin(np.radians(beta))])
        s=solve_ivp(rhs,(1.0,T[2]),[Dt(1.0),0.0],args=(n,),events=ev,rtol=1e-10,atol=1e-14,max_step=2e-3)
        te=s.t_events[0][0]; xe=s.y_events[0][0]; qe=q_of_r(Rr,te); assert abs(r_of(qe,te)-Rr)<1e-6*Rr+1e-9
        re,dre=lay(qe,te); rn,drn=lay(qe,1.0); S=((rn**2*drn)/(re**2*dre))**(1/3)
        u=xe/Rr; Ve=v_of_r(rn,1.0)*u; Vu=np.array([v_of_r(Dt(1.0),1.0),0]); vr=(Ve-Vu)/C
        los=rn*u-[Dt(1.0),0]; los/=np.hypot(*los); br=vr@los; b2=vr@vr
        out.append(dict(beta=beta,te=te,S=S,dop=1+br,gam=1/np.sqrt(1-b2) if b2<1 else np.nan,b=np.sqrt(b2)))
    return out
meanT=lambda Rr: np.mean([3000/r['S'] for r in rays(Rr,30)])-2.725
lo,hi=0.005/DG,0.1/DG

RR=brentq(lambda x: meanT(np.exp(x)),np.log(lo),np.log(hi),xtol=1e-4)
RR=np.exp(RR); print('R_rec Mpc',RR*DG*1e3,flush=True)
def st(Tv,b):
    A=np.c_[np.ones_like(b),np.cos(np.radians(b)),np.sin(np.radians(b))]; c,*_=np.linalg.lstsq(A,Tv,rcond=None); r=Tv-A@c
    return dict(mean=float(Tv.mean()),dipole_K=float(np.hypot(c[1],c[2])),left_rms_K=float(r.std()),left_1e5=float(r.std()/Tv.mean()/1e-5))
res=dict(R_rec_Mpc=RR*DG*1e3,trec_Gyr=float(np.interp(RR,XR[:,i1],T)*t0),obs_shift_vs_code_taui0p01_Gpc=float(DIFF[0]*DG),
         phi_c2_Rrec=phi(RR)/C**2,phi_c2_obs=phi(Dt(1.0))/C**2)
for step in (1,6):
    R=rays(RR,step); g=lambda k: np.array([r[k] for r in R]); b=g('beta'); S,d,gm=g('S'),g('dop'),g('gam')
    res[f'{step}deg']=dict(terms={k:[float(v.min()),float(v.mean()),float(v.max())] for k,v in dict(stretch=S,doppler=d,gamma=gm,beta=g('b')).items()},
        stretch_only=st(3000/S,b),stretch_dop=st(3000/(S*d),b),stretch_dop_gam=st(3000/(S*d*gm),b),
        share=dict(stretch=float(np.mean(np.log(S))/np.mean(np.log(S*d*gm))),motion=float(np.mean(np.log(d*gm))/np.mean(np.log(S*d*gm))),grav=0.0))
    if step==1: np.save('cmbann_keep.npy',np.array([b,S,d,gm,g('te')]))
# annulus: us at centre, bounding layers D±R now (R=ct*), traced back
tr=res['trec_Gyr']/t0; Rn=C*(1-tr); qi,qo=q_of_r(Dt(1)-Rn,1.0),q_of_r(Dt(1)+Rn,1.0)
assert np.isclose(r_of(qi,1)+r_of(qo,1),2*Dt(1)) and np.isclose(r_of(qo,1)-r_of(qi,1),2*Rn)
ann=dict(q_in=qi,q_out=qo,now=[r_of(qi,1)*DG,r_of(qo,1)*DG],at_trec=[r_of(qi,tr)*DG,r_of(qo,tr)*DG],D_trec=Dt(tr)*DG)
V=lambda a,b_: b_**3-a**3; ann['V_ratio']=V(*ann['now'])/V(*ann['at_trec']); ann['stretch']=ann['V_ratio']**(1/3); res['annulus']=ann
json.dump(res,open('cmbann.json','w'),indent=1,default=float); print(json.dumps(res,indent=1,default=float))
