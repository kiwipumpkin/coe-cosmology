import sys,json,pickle,numpy as np,multiprocessing as mp;sys.path.insert(0,'..');from rings import Flow
from scipy.optimize import brentq
import rings
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline
class FlowW(Flow):
    def __init__(s,k,GM,tmax=6.0):
        s.k,s.GM=k,GM;s.r0=np.geomspace(2.05,60,500);s.t=np.linspace(0,tmax,601);Rr=[];Vv=[]
        for r0 in s.r0:
            so=solve_ivp(lambda t,y:[y[1],-GM/y[0]**2],(0,tmax),[r0,k*r0],t_eval=s.t,rtol=1e-11,atol=1e-12);Rr.append(so.y[0]);Vv.append(so.y[1])
        s.R=RectBivariateSpline(s.r0,s.t,np.array(Rr),kx=5,ky=5);s.V=RectBivariateSpline(s.r0,s.t,np.array(Vv),kx=5,ky=5)
def grad(r0,kind,q):   # initial temperature profile T0(r0)/T0(3)
    if kind=='pow':return (r0/3.0)**(-q)
    if kind=='exp':return np.exp(-(r0-3.0)/q)
def sky_sr(F,c,tobs,trec,kind,q,N=180):
    # emission when T=3000: T = T0(r0)/aloc => aloc(r0,t)/grad(r0) = aloc(3,trec)   (observer's own element recombines at trec)
    X0=F.aloc(3.0,trec);xo=np.array([F.R(3.0,tobs,grid=False),0.]);v=F.V(3.0,tobs,grid=False)/c;g=1/np.sqrt(1-v*v);res=[]
    for th in np.radians(np.arange(N)*360/N):
        n=np.array([np.cos(th),np.sin(th)])            # observer-frame sight direction
        def ev(tau):  # observer-frame lookback (time units); lab event
            tp=-tau;xp=n*c*tau;t=tobs+g*(tp+v*xp[0]/c);x=xo+np.array([g*(xp[0]+v*c*tp),xp[1]]);return t,x
        def h(tau):
            t,x=ev(tau)
            if t<0:return np.nan
            r=np.hypot(*x);q0=F.r0_of(r,t)
            if not np.isfinite(q0):return np.nan   # medium extends beyond the plotted disc (all r0 in grid)
            return F.aloc(q0,t)/grad(q0,kind,q)-X0
        taus=np.linspace(0,tobs*1.2,1501);hv=np.array([h(s) for s in taus]);ok=np.isfinite(hv)
        idx=np.where(ok[:-1]&ok[1:]&(np.sign(hv[:-1])!=np.sign(hv[1:])))[0]
        if len(idx)==0:res.append((np.nan,)*4);continue
        i=idx[0];tau=brentq(h,taus[i],taus[i+1]);t,x=ev(tau);r=np.hypot(*x);q0=F.r0_of(r,t)
        u=F.V(q0,t,grid=False)/c*x/r;d=1-u[0]*v;up=np.array([(u[0]-v)/d,u[1]/(g*d)]);ge=1/np.sqrt(1-up@up)
        res.append((1/(ge*(1-up@(-n))),t,tau,q0))
    return np.array(res)
CASES=[(1.6,12.0),(0.8,2.0),(1.6,0.0)];GR=[('pow',0.0),('pow',0.5),('pow',1.0),('pow',2.0),('exp',1.0)];C=20.0;TOBS=3.0;TREC=2.95
def job(a):
    (k,G),(kind,q)=a;F=FlowW(k,G);return a,sky_sr(F,C,TOBS,TREC,kind,q)
if __name__=='__main__':
    jobs=[(cs,gr) for cs in CASES for gr in GR]
    with mp.Pool(8) as P:R=dict(P.map(job,jobs))
    pickle.dump(R,open('v2.pkl','wb'));rows=[]
    for ((k,G),(kind,q)),x in R.items():
        T=x[:,0];ok=np.isfinite(T);d=T/np.nanmean(T)-1;z=np.nan_to_num(d);m=np.fft.rfft(z)/len(z)
        rows.append(dict(k=k,GM=G,grad=kind,q=q,missing=int((~ok).sum()),rms=float(np.sqrt(np.nanmean(d**2))),range=float(np.nanmax(d)-np.nanmin(d)),dipole=float(2*abs(m[1])),
          dipole_sign=float(np.sign(m[1].real)),quad=float(2*abs(m[2])),te_spread=float(np.nanmax(x[:,1])-np.nanmin(x[:,1])),te_0=float(x[0,1]),te_90=float(x[45,1]),te_180=float(x[90,1])))
        print(rows[-1],flush=True)
    json.dump(rows,open('v2_stats.json','w'),indent=1)
