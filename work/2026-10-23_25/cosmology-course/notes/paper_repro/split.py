import numpy as np,json
exec(open('cmbann.py').read().split("meanT=lambda")[0])
from scipy.optimize import brentq
d=json.load(open('cmbann.json')); b,S1,dop,gm,te=np.load('cmbann_keep.npy')
def geom(Rr):
    tr=np.interp(Rr,XR[:,i1],T); Rn=C*(1-tr); D1=Dt(1.0)
    qi,qo=q_of_r(D1-Rn,1.0),q_of_r(D1+Rn,1.0)
    V=lambda a,c,t: r_of(c,t)**3-r_of(a,t)**3
    si=(V(qi,1.0,1.0)/V(qi,1.0,tr))**(1/3); so=(V(1.0,qo,1.0)/V(1.0,qo,tr))**(1/3); sa=(V(qi,qo,1.0)/V(qi,qo,tr))**(1/3)
    return tr,qi,qo,si,so,sa
cb=np.cos(np.radians(b))
def Tsplit(Rr):
    tr,qi,qo,si,so,sa=geom(Rr); S=si*(1-cb)/2+so*(1+cb)/2; return 3000/S,S
RR=np.exp(brentq(lambda x: Tsplit(np.exp(x))[0].mean()-2.725,np.log(0.005/DG),np.log(0.1/DG),xtol=1e-6))
tr,qi,qo,si,so,sa=geom(RR); Ts,Ss=Tsplit(RR)
def st(Tv):
    A=np.c_[np.ones_like(b),cb,np.sin(np.radians(b))]; c,*_=np.linalg.lstsq(A,Tv,rcond=None); r=Tv-A@c
    return dict(mean=float(Tv.mean()),dipole_K=float(np.hypot(c[1],c[2])),left_rms_K=float(r.std()),left_1e5=float(r.std()/Tv.mean()/1e-5))
b6=slice(None,None,6)
out=dict(R_rec_Mpc=RR*DG*1e3,trec_Myr=tr*t0*1e3,q_in=qi,q_out=qo,stretch_inner=si,stretch_outer=so,stretch_whole=sa,
  split_1deg=st(Ts),split_6deg=None,single_ray_1deg=d['1deg']['stretch_only'],
  subshells_now_Gpc=[r_of(qi,1)*DG,Dt(1)*DG,r_of(qo,1)*DG],subshells_trec_Mpc=[r_of(qi,tr)*DG*1e3,Dt(tr)*DG*1e3,r_of(qo,tr)*DG*1e3])
bb=b[b6];T6=Ts[b6];A=np.c_[np.ones_like(bb),np.cos(np.radians(bb)),np.sin(np.radians(bb))];c,*_=np.linalg.lstsq(A,T6,rcond=None);r=T6-A@c
out['split_6deg']=dict(mean=float(T6.mean()),dipole_K=float(np.hypot(c[1],c[2])),left_rms_K=float(r.std()),left_1e5=float(r.std()/T6.mean()/1e-5))
assert r_of(qi,1)<Dt(1)<r_of(qo,1) and np.isclose(r_of(qo,1)-r_of(qi,1),2*C*(1-tr),rtol=1e-3)
json.dump(out,open('split.json','w'),indent=1,default=float); np.save('split_T.npy',np.array([b,Ts,3000/S1])); print(json.dumps(out,indent=1,default=float))
