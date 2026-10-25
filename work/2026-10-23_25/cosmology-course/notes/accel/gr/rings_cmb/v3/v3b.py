# Recombination sphere centred on CoE at r_rec < 2 (inside shell 10). Elements crossed it outward at t<0 (flow integrated backward).
import sys,json,pickle,numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq
R0=np.linspace(2.0,4.0,201)
class FlowB:
    def __init__(s,k,GM):
        s.k,s.GM=k,GM;s.sol={}
        for r0 in R0:
            f=lambda t,y:[y[1],-GM/y[0]**2];ev=lambda t,y:y[0]-0.02;ev.terminal=True
            b=solve_ivp(f,(0,-10),[r0,k*r0],events=ev,dense_output=True,rtol=1e-11,atol=1e-12)
            a=solve_ivp(f,(0,8),[r0,k*r0],dense_output=True,rtol=1e-11,atol=1e-12);s.sol[r0]=(b,a)
    def state(s,r0,t):
        b,a=s.sol[r0];so=a if t>=0 else b
        if t<so.t[-1] if t<0 else t>so.t[-1]:return None
        return so.sol(t)
    def tcross(s,r0,R):
        b,a=s.sol[r0];f=lambda t:b.sol(t)[0]-R
        return brentq(f,b.t[-1],0) if f(b.t[-1])<0<f(0) else np.nan
def sky(F,c,tobs,rrec,N=180):
    # crossing events: element r0 at |x|=rrec at time tc(r0); emission event lies on circle |x|=rrec at time t in [min tc, max tc]
    tc=np.array([F.tcross(r0,rrec) for r0 in R0]);ok=np.isfinite(tc)
    r0_of_t=interp1d(tc[ok],R0[ok],bounds_error=False)   # monotonic: which element is at rrec at time t
    so=F.state(3.0,tobs);xo=np.array([so[0],0.]);v=so[1]/c;g=1/np.sqrt(1-v*v);out=[]
    for th in np.radians(np.arange(N)*360/N):
        n=np.array([np.cos(th),np.sin(th)])
        def ev(tau):tp=-tau;xp=n*c*tau;return tobs+g*(tp+v*xp[0]/c),xo+np.array([g*(xp[0]+v*c*tp),xp[1]])
        h=lambda tau:np.hypot(*ev(tau)[1])-rrec
        taus=np.linspace(0,15,6001);hv=np.array([h(s) for s in taus]);idx=np.where(np.sign(hv[:-1])!=np.sign(hv[1:]))[0];got=None
        for i in idx:
            tau=brentq(h,taus[i],taus[i+1]);t,x=ev(tau);r0=float(r0_of_t(t))
            if np.isfinite(r0):got=(t,x,r0);break
        if got is None:out.append((np.nan,)*3);continue
        t,x,r0=got;st=F.sol[R0[np.argmin(abs(R0-r0))]][0].sol(t);u=st[1]/c*x/rrec
        if u@u>=1:out.append((np.nan,)*3);continue
        d=1-u[0]*v;up=np.array([(u[0]-v)/d,u[1]/(g*d)]);ge=1/np.sqrt(1-up@up);out.append((1/(ge*(1-up@(-n))),t,r0))
    return np.array(out),tc
if __name__=='__main__':
    rows=[];R={}
    for k,G in [(1.6,12.0),(1.6,0.0)]:
        F=FlowB(k,G)
        for rrec in (0.5,1.0,1.5,1.9):
            for c in (20.0,50.0):
                for tobs in (0.0,1.0,3.0):
                    x,tc=sky(F,c,tobs,rrec);T=x[:,0];ok=np.isfinite(T);row=dict(k=k,GM=G,r_rec=rrec,c=c,tobs=tobs,valid=int(ok.sum()),
                        cross_shell1_outer=float(tc[-1]),cross_shell10_inner=float(tc[0]),order_ok=bool(np.all(np.diff(tc[np.isfinite(tc)])>0)))
                    if ok.sum()>=3:
                        d=T/np.nanmean(T)-1;z=np.nan_to_num(d);m=np.fft.rfft(z)/len(z)
                        row.update(rms=float(np.sqrt(np.nanmean(d**2))),range=float(np.nanmax(d)-np.nanmin(d)),dipole=float(2*abs(m[1])),dipole_sign=float(np.sign(m[1].real)),quad=float(2*abs(m[2])),
                          te_min=float(np.nanmin(x[:,1])),te_max=float(np.nanmax(x[:,1])),te_0=float(x[0,1]),te_90=float(x[45,1]),te_180=float(x[90,1]))
                        R[(k,G,rrec,c,tobs)]=x
                    rows.append(row);print({a:(round(b,5) if isinstance(b,float) else b) for a,b in row.items()},flush=True)
    json.dump(rows,open('v3b_stats.json','w'),indent=1);pickle.dump(R,open('v3b.pkl','wb'))
