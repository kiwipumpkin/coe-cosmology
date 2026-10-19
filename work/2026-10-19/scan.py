import sys,json,numpy as np,v11,base as V
from scipy.optimize import brentq
M=0.008532635720463873;R=float(sys.argv[1]);F=v11.field(M,R);wo,t0,T,W,S=F
ps=np.arange(0,181,5.0);o,_=v11.trace(F,R,ps);ok=np.isfinite(o[:,1])
rt=np.array([V.Rfun(t0,w,M) for w in o[:,4]]);phi=o[:,5];xo=V.Rfun(t0,wo,M)
tc=brentq(lambda t:V.Rfun(t,wo,M)-R,T[0]*1.01,t0);rv=t0-tc
d=np.hypot(rt*np.cos(phi)-xo,rt*np.sin(phi));d[~ok]=np.nan;i=np.nanargmin(d)
G=4.2/rv
json.dump(dict(miss=ps[~ok].tolist(),d180=float(d[-1]) if ok[-1] else None,R=R,ok=int(ok.sum()),n=len(ps),circle=rv,minred=float(d[i]),ang=float(ps[i]),D=G,age=t0*G*3.26156,tc=tc,t0=t0),open(f'scan_{R}.json','w'))
