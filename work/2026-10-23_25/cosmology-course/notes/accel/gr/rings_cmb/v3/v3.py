import sys,json,pickle,numpy as np;sys.path.insert(0,'..');sys.path.insert(0,'../v2')
from v2 import FlowW
from scipy.optimize import brentq
RREC=float(sys.argv[1]) if len(sys.argv)>1 else 4.0
def tcross(F,r0,R):   # time element r0 reaches radius R
    f=lambda t:F.R(r0,t,grid=False)-R
    return brentq(f,0,6) if f(0)<=0<=f(6) else np.nan
def sky(F,c,tobs,N=180):
    xo=np.array([F.R(3.0,tobs,grid=False),0.]);v=F.V(3.0,tobs,grid=False)/c;g=1/np.sqrt(1-v*v);out=[]
    for th in np.radians(np.arange(N)*360/N):
        n=np.array([np.cos(th),np.sin(th)])
        def ev(tau):tp=-tau;xp=n*c*tau;return tobs+g*(tp+v*xp[0]/c),xo+np.array([g*(xp[0]+v*c*tp),xp[1]])
        h=lambda tau:np.hypot(*ev(tau)[1])-RREC
        taus=np.linspace(0,3,3001);hv=np.array([h(s) for s in taus]);idx=np.where(np.sign(hv[:-1])!=np.sign(hv[1:]))[0]
        got=None
        for i in idx:
            tau=brentq(h,taus[i],taus[i+1]);t,x=ev(tau)
            if t<0:break
            r0=F.r0_of(RREC,t)
            if np.isfinite(r0) and 2.0<=r0<=4.0:got=(tau,t,x,r0);break
        if got is None:out.append((np.nan,)*3);continue
        tau,t,x,r0=got;u=F.V(r0,t,grid=False)/c*x/RREC;d=1-u[0]*v;up=np.array([(u[0]-v)/d,u[1]/(g*d)]);ge=1/np.sqrt(1-up@up)
        out.append((1/(ge*(1-up@(-n))),t,r0))
    return np.array(out)
if __name__=='__main__':
    rows=[];R={}
    for k,G in [(1.6,12.0),(1.6,0.0)]:
        F=FlowW(k,G)
        ts=[tcross(F,r0,RREC) for r0 in np.arange(2.0,4.01,0.2)]   # shell edges (outer=4 -> shell 1 ... inner=2 -> shell 10)
        mono=bool(np.all(np.diff(np.array(ts)[::-1])>=0)) if True else None
        # no-crossing check: R(r0,t) increasing in r0 for all t
        tt=np.linspace(0,6,61);rr=np.linspace(2.05,4.0,40);Rg=np.array([[F.R(r,t,grid=False) for r in rr] for t in tt]);nocross=bool(np.all(np.diff(Rg,axis=1)>0))
        tc3=tcross(F,3.0,RREC)
        print(k,G,'edge crossing times r0=2.0..4.0:',np.round(ts,4),'nocross',nocross,'observer crosses at',tc3)
        for c in (20.0,50.0):
            for tobs in ([0.9*tc3,1.0,2.0] if np.isfinite(tc3) else [1.0,2.0]):
                x=sky(F,c,tobs);T=x[:,0];ok=np.isfinite(T)
                if ok.sum()<3:rows.append(dict(k=k,GM=G,c=c,tobs=tobs,valid=int(ok.sum())));print(rows[-1]);continue
                d=T/np.nanmean(T)-1;z=np.nan_to_num(d);m=np.fft.rfft(z)/len(z)
                rows.append(dict(k=k,GM=G,c=c,tobs=float(tobs),valid=int(ok.sum()),valid_angles=[float(a) for a in np.arange(0,360,2.0)[ok][[0,-1]]],rms=float(np.sqrt(np.nanmean(d**2))),range=float(np.nanmax(d)-np.nanmin(d)),
                  dipole_ring=float(2*abs(m[1])),dipole_sign=float(np.sign(m[1].real)),quad_ring=float(2*abs(m[2])),te_min=float(np.nanmin(x[:,1])),te_max=float(np.nanmax(x[:,1])),
                  shells_seen=sorted(set(int(np.clip(10-np.floor((r-2.0)/0.2),1,10)) for r in x[ok,2]))))
                R[(k,G,c,float(tobs))]=x;print(rows[-1],flush=True)
        rows.append(dict(k=k,GM=G,edge_cross_times_r0_2_to_4=[float(a) for a in ts],no_shell_crossing=nocross,observer_crosses_Rrec_at=float(tc3)))
    json.dump(rows,open(f'v3_stats_R{RREC}.json','w'),indent=1);pickle.dump(R,open(f'v3_R{RREC}.pkl','wb'))
