# Job 1: shell built as an exact boost of a round shell at rest (radius 1, expanding at u in its own frame).
# All particles start on the moving frame's slice tau = v0 (so every start event has t >= 0); lab position/time:
#   t_s = g v0 (1+e_x),  X_s = (g (e_x + v0^2), e_y);  observer: t_s = g v0, X_s = (g v0^2, 0)
import sys,pickle; sys.path.insert(0,'.'); from grcore import *
from scipy.integrate import solve_ivp
import multiprocessing as mp
def starts(v0):
    g=1/np.sqrt(1-v0*v0); E=[]
    for t in th:
        e=np.array([np.cos(t),np.sin(t)]); E.append((g*v0*(1+e[0]), np.array([g*(e[0]+v0*v0),e[1]]), e))
    return (g*v0, np.array([g*v0*v0,0.0])), E
def geod2(M,X0,w0,t0,t1):
    at,h=geodesic(M,X0,w0,t0,t1)
    return (lambda t: None if t<t0-1e-12 else at(t)), h
def observe2(M,Xo,Uo,ems,to):
    ro,po=polar(Xo); out=[]; prev=None
    for em,ts in ems:
        g=lambda s: np.linalg.norm(Xo-em(s)[0])-(to-s)
        if g(ts)>0: out.append(np.nan); continue   # light from this emitter's start has not arrived yet
        hi=to-1e-6
        while em(hi) is None: hi=ts+(hi-ts)*0.999
        s0=brentq(g,ts,hi) if g(hi)<0 is False or g(hi)>0 else hi
        Xe=em(s0)[0]; d=Xe-Xo; dr=d@np.array([np.cos(po),np.sin(po)]); dp=d@np.array([-np.sin(po),np.cos(po)])
        fo=1-2*M/ro; b0=np.arctan2(-dp,-dr/np.sqrt(fo))
        def res(z):
            b,s=z
            if not (ts<=s<=to): return [10+abs(s),10]
            so,L=ray(M,Xo,b,to,ts-1)
            if so.t[-1]>s: return [5,5]
            e=em(s)
            if e is None: return [7,7]
            r,p,_=so.sol(s); return cart(r,p)-e[0]
        sol=None
        for guess in ([b0,s0],prev):
            if guess is None: continue
            z=fsolve(res,guess,xtol=1e-12)
            if np.max(np.abs(res(z)))<1e-8: sol=z;break
        if sol is None: out.append(np.nan); continue
        prev=sol; b,s=sol; so,L=ray(M,Xo,b,to,ts-1); r,p,pr=so.sol(s)
        out.append(kdotu(M,ro,np.cos(b),L,Uo)/kdotu(M,r,pr,L,em(s)[1]))
    return np.array(out)
def gr_boost(mu,u,v0,to):
    (ots,oX),E=starts(v0); ems=[];hz=0
    for ts,X,e in E:
        at,h=geod2(mu,X,relsum(v0,u,e),ts,to); ems.append((at,ts)); hz+=h
    oat,_=geod2(mu,oX,np.array([v0,0]),ots,to); Xo,Uo=oat(to)
    return observe2(mu,Xo,Uo,ems,to),hz
def weak_boost(mu,u,v0,to):
    def tr(p0,w,t0):
        def f(t,y):
            p=y[2:]; v=p/np.sqrt(1+p@p); F=-mu*(y[:2]-Mp)/np.linalg.norm(y[:2]-Mp)**3
            return [*v,*F]
        g=1/np.sqrt(1-w@w); return solve_ivp(f,(t0,to),[*p0,*(g*w)],dense_output=True,rtol=1e-11,atol=1e-13).sol
    (ots,oX),E=starts(v0); yo=tr(oX,np.array([v0,0]),ots)(to); xo=yo[:2]; vo=yo[2:]/np.sqrt(1+yo[2:]@yo[2:]); out=[]
    for ts,X,e in E:
        S=tr(X,relsum(v0,u,e),ts); P=lambda s:S(s)[:2]
        g=lambda s: np.linalg.norm(xo-P(s))-(to-s)
        if g(ts)>0: out.append(np.nan); continue
        te=brentq(g,ts,to-1e-6); y=S(te); xe=y[:2]; ve=y[2:]/np.sqrt(1+y[2:]@y[2:]); n=(xo-xe)/np.linalg.norm(xo-xe)
        re=np.linalg.norm(xe-Mp); ro=np.linalg.norm(xo-Mp)
        if mu and re<=2*mu: out.append(np.nan); continue
        gg=np.sqrt((1-2*mu/re)/(1-2*mu/ro)) if mu else 1.0
        out.append(1/np.sqrt(1-vo@vo)*(1-vo@n)/(1/np.sqrt(1-ve@ve)*(1-ve@n))*gg)
    return np.array(out)
def job(c):
    to,mu,u,v0=c; return c,weak_boost(mu,u,v0,to),gr_boost(mu,u,v0,to)
if __name__=='__main__':
    C=[(to,mu,u,v0) for v0,to in ((0.1,2.3),(0.3,2.3),(0.5,2.3),(0.9,11.0),(0.99,106.0)) for u in (0,0.9999) for mu in (0.0,0.1)]
    with mp.Pool(8) as P: R=dict((c,(w,g)) for c,w,g in P.map(job,C))
    pickle.dump(R,open('boost.pkl','wb'))
    s=lambda r:("%.4f %.4f %.4f"%(stats(norm(r))[0],stats(norm(r))[1],stats(norm(r))[3])) if not np.isnan(r).any() else "missing %d"%np.isnan(r).sum()
    for c,(w,g) in sorted(R.items()): print(c,"weak",s(w),"| GR",s(g[0]),"hz",g[1],flush=True)
