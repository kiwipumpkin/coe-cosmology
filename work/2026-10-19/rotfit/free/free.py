import sys;sys.path.insert(0,'..')
import numpy as np,healpy as hp,json
from scipy.optimize import minimize
from fit import sky,vec,lb,ang,dip,quad,V,T0,T1,psi,PH,PQ,PR,nside
PC=-PH
def hang(a,b):o=ang(a,b);return min(o,180-o)
# fixed no-motion axial dipole (monopole + cos fit on the no-motion profile, as in pec.py)
p=np.radians(psi);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(.5))/4;w/=w.sum()
d0=T1(psi)/np.sum(w*T1(psi))-1;X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w)
c0,c1=np.linalg.lstsq(X*sw[:,None],d0*sw,rcond=None)[0]
def resid(a,v,n=V):
    m=sky(a,v);return m-c0-c1*(n@a)
def q4(x):return vec(x[0],x[1]),vec(x[2],x[3])   # anti-CoE axis, motion
def costA(x):a,v=q4(x);return ang(dip(sky(a,v)),PH)
def costB(x):a,v=q4(x);r=resid(a,v);return ang(dip(r),PH)+hang(quad(r)[0],PQ)
def mfit(cost,starts):
    best=None
    for s in starts:
        r=minimize(cost,s,method='Nelder-Mead',options=dict(xatol=.05,fatol=1e-4,maxiter=3000))
        if best is None or r.fun<best.fun:best=r
    return best
out={'roll':'null: model sky depends only on angle from the axis (quadrupole eigenvalues equal to 3e-4 relative), so roll is not a knob'}
lPC=lb(PC);lPQ=lb(-PQ)
rA=mfit(costA,[[*lPC,*lb(PH)],[*lPC,254.7,46.2],[0,15,254.7,46.2]])
stB=[]
for l in (0,90,180,270):
    for qs in (lb(PQ),lb(-PQ)): stB.append([*qs,264,48]);stB.append([*qs,l,0])
rB=mfit(costB,stB)
def report(x,kind):
    a,v=q4(x);m=sky(a,v) if kind=='A' else resid(a,v);d=dip(m);q,e=quad(m)
    al=lb(a);cl=lb(-a)
    return dict(antiCoE_lb=al,CoE_lb=cl,motion_lb=lb(v),dipole_hot_lb=lb(d),hot_offset_deg=float(ang(d,PH)),cold_offset_deg=float(ang(-d,PC)),
      quad_axis_lb=lb(q),quad_offset_deg=float(hang(q,PQ)),dipole_amp_uK=float(np.linalg.norm(d)*T0),dipole_amp_over_Planck=float(np.linalg.norm(d)*T0/3362),
      quad_rms_uK=float(np.sqrt(np.mean(hp.alm2map(np.where(np.arange(6)>=3,hp.map2alm(m,lmax=2),0),nside)**2))*T0) if False else None,
      range_over_Planck_ring=float(np.ptp(m)/PR),x=list(map(float,x)))
out['A_full_dipole']=report(rA.x,'A');out['B_axial_dipole_removed']=report(rB.x,'B')
# quadrupole rms per case
for k,kind in (('A_full_dipole','A'),('B_axial_dipole_removed','B')):
    a,v=q4(out[k]['x']);m=sky(a,v) if kind=='A' else resid(a,v)
    al=hp.map2alm(m,lmax=2,iter=3);ll=hp.Alm.getlm(2)[0];al[ll!=2]=0;out[k]['quad_rms_uK']=float(np.std(hp.alm2map(al,nside))*T0)
print(json.dumps(out,indent=1,default=float))
json.dump(out,open('free_fit.json','w'),indent=1,default=float)
