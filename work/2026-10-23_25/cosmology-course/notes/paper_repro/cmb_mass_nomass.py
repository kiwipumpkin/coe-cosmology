import numpy as np,json,sys
from scipy.optimize import brentq
MUSCALE=float(sys.argv[1])  # 1 = fitted central mass, 0 = no mass (coasting)
s=open('cmbann.py').read().split("meanT=lambda")[0]
s=s.replace("exec(src)","src=src.replace(\"mu,MT,RG,S=rr['mu'],rr['MT'],pc.RG,rr['S']\",\"mu,MT,RG,S=rr['mu']*%r,rr['MT'],pc.RG,rr['S']\"); exec(src)"%MUSCALE)
exec(s)
b=np.arange(360.); cb=np.cos(np.radians(b))
def geom(Rr):
    tr=np.interp(Rr,XR[:,i1],T); Rn=C*(1-tr); D1=Dt(1.0)
    qi,qo=q_of_r(D1-Rn,1.0),q_of_r(D1+Rn,1.0)
    V=lambda a,c,t: r_of(c,t)**3-r_of(a,t)**3
    si=(V(qi,1.0,1.0)/V(qi,1.0,tr))**(1/3); so=(V(1.0,qo,1.0)/V(1.0,qo,tr))**(1/3)
    return tr,si,so
def Tsplit(Rr):
    tr,si,so=geom(Rr); return 3000/(si*(1-cb)/2+so*(1+cb)/2)
RR=np.exp(brentq(lambda x: Tsplit(np.exp(x)).mean()-2.725,np.log(0.002/DG),np.log(0.2/DG),xtol=1e-8))
Ts=Tsplit(RR); tr,si,so=geom(RR); dT=(Ts-Ts.mean())/Ts.mean()
A=np.c_[np.ones_like(b),cb,np.sin(np.radians(b))]; c,*_=np.linalg.lstsq(A,dT,rcond=None)
assert abs(Ts.mean()-2.725)<1e-6
out=dict(mu=float(mu),R_rec_Mpc=RR*DG*1e3,D_now_Gpc=Dt(1.0)*DG,stretch_inner=si,stretch_outer=so,
 full_amp_1e5=float((dT.max()-dT.min())/2/1e-5),dipole_1e5=float(np.hypot(c[1],c[2])/1e-5),
 hot_beta=float(np.degrees(np.arctan2(c[2],c[1]))%360),hot_l=float((np.degrees(np.arctan2(c[2],c[1]))+15.2)%360))
np.save('cmb_mass%g.npy'%MUSCALE,np.array([b,dT])); json.dump(out,open('cmb_mass%g.json'%MUSCALE,'w'),indent=1); print(out)
