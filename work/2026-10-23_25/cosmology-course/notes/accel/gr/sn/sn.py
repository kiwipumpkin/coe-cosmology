# Luminosity distance in the LTB valley model, comoving observer at the model's Earth position (D,0); rays in the plane through the CoE axis.
import sys,json,pickle,os,warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/valley'); sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/ltb')
import numpy as np, toy
from ltbcore import LTB, photon, edge
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline, CubicSpline
RV=14.30048287260589
class LTBE(LTB):
    def __init__(s,Mfun,D,Rdot,rmin,rmax=60.0,nr=500):
        s.D=D; M=Mfun; s.Mf=M; s.Mc=0.0
        ED=(Rdot**2-2*M(D)/D)/2; s.t1=s.age(D,M(D),ED)
        rg=np.linspace(rmin,rmax,nr); Eg=np.array([s.solveE(r,M(r)) for r in rg]); s.Es=CubicSpline(rg,Eg)
        tg=np.unique(np.r_[np.geomspace(0.03*s.t1,s.t1,400), s.t1+np.geomspace(1e-3,20,60)])
        Rg=np.empty((nr,len(tg)))
        for i,(r,E) in enumerate(zip(rg,Eg)):
            Mr=M(r); f=lambda t,R:[np.sqrt(max(2*Mr/R[0]+2*E,1e-30))]
            so=solve_ivp(f,(s.t1,tg[-1]),[r],t_eval=tg[tg>=s.t1],rtol=1e-11,atol=1e-12)
            sb=solve_ivp(f,(s.t1,tg[0]),[r],t_eval=tg[tg<s.t1][::-1],rtol=1e-11,atol=1e-12)
            Rg[i]=np.r_[sb.y[0][::-1],so.y[0]]
        s.Rs=RectBivariateSpline(rg,tg,Rg,kx=5,ky=5); s.rg=rg; s.tg=tg; s.Eg=Eg
        s.Rr=s.Rs.partial_derivative(1,0); s.Rt_=s.Rs.partial_derivative(0,1); s.Rrt=s.Rs.partial_derivative(1,1); s.Rrr=s.Rs.partial_derivative(2,0)
        s.minRp=float(s.Rs(rg,tg,dx=1).min())
def build(i,rmin):
    p=json.load(open('/workspace/cosmology-course/notes/accel/gr/valley/valley_points.json'))[i]; fac,Dt=toy.make(p)
    return LTBE(fac.args[0],Dt,toy.HR*Dt,rmin),Dt,p
def ray(L,D,beta,dB=1e-4):
    tl=0.035*L.t1; out=[]
    for b in (beta,beta+dB):
        so=photon(L,D,0.0,b,L.t1,tl); out.append(so)
    t=np.linspace(L.t1,max(out[0].t[-1],out[1].t[-1]),3000)[1:]
    y=out[0].sol(t); y2=out[1].sol(t)
    R,Rp,Rt,G,_,_=L.fields(t,y[0]) if False else (None,)*6
    res=[]
    for k in range(len(t)):
        r,ps,vr,vp,lnk=y[:,k]; R,Rp,Rt,G,_,_=L.fields(t[k],r)
        xi=np.array([y2[0,k]-r,(y2[1,k]-ps)])/dB         # deviation at fixed coordinate time
        g=np.diag([G,R*R]); v=np.array([vr,vp]); vn=v/np.sqrt(v@g@v)
        xp=xi-(xi@g@vn)*vn; Dpar=np.sqrt(xp@g@xp)
        Dperp=R*abs(np.sin(ps))/abs(np.sin(beta))
        res.append((np.exp(lnk)-1,np.sqrt(Dpar*Dperp),t[k],r))
    a=np.array(res); a=a[np.argsort(a[:,0])]; return a  # z, d_A, t, r
if __name__=='__main__':
    i=int(sys.argv[1]); rmin=float(sys.argv[2]); L,Dt,p=build(i,rmin)
    betas=np.radians([1.0]+list(range(15,180,15))+[179.0]); out=dict(i=i,Dt=Dt,t1=L.t1,minRp=L.minRp,rays={})
    print('built',i,L.t1,L.minRp,flush=True)
    for b in betas:
        a=ray(L,Dt,b); out['rays'][float(np.degrees(b))]=a; pickle.dump(out,open(f'dl_{i}.pkl','wb'))
        print(np.degrees(b),'zmax %.3f'%a[:,0].max(),flush=True)
