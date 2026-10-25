import sys,pickle; sys.path.insert(0,'.'); from grcore import *
import multiprocessing as mp
def weak_T(xo,vo,emit,to,mu,lo):
    out=[]
    for t in th:
        e=np.array([np.cos(t),np.sin(t)])
        try: te=brentq(lambda s: np.linalg.norm(xo-emit(s,e)[0])-(to-s),lo,to-1e-6)
        except Exception: return None
        xe,ve=emit(te,e); n=(xo-xe)/np.linalg.norm(xo-xe)
        re=np.linalg.norm(xe-Mp); ro=np.linalg.norm(xo-Mp)
        if mu and (re<=2*mu or ro<=2*mu) : return None
        ge,go=1/np.sqrt(1-ve@ve),1/np.sqrt(1-vo@vo)
        g=np.sqrt((1-2*mu/re)/(1-2*mu/ro)) if mu else 1.0
        out.append(go*(1-vo@n)/(ge*(1-ve@n))*g)
    return norm(np.array(out))
def weak_coast(mu,u,v0,to=3.0):
    def tr(p0,w):
        def f(t,y):
            p=y[2:]; v=p/np.sqrt(1+p@p); F=-mu*(y[:2]-Mp)/np.linalg.norm(y[:2]-Mp)**3
            return [*v,*F]
        g=1/np.sqrt(1-w@w); return solve_ivp(f,(0,to),[*p0,*(g*w)],dense_output=True,rtol=1e-11,atol=1e-13).sol
    yo=tr([0,0],np.array([v0,0]))(to); vo=yo[2:]/np.sqrt(1+yo[2:]@yo[2:]); cache={}
    def emit(s,e):
        key=(round(e[0],12),round(e[1],12))
        if key not in cache: cache[key]=tr(list(e),relsum(v0,u,e))
        y=cache[key](s); return y[:2],y[2:]/np.sqrt(1+y[2:]@y[2:])
    return weak_T(yo[:2],vo,emit,to,mu,0.0)
def boostpos(v,u,e):
    G=1/np.sqrt(1-v*v); tau=(0/G - v*e[0])/(1+v*u*e[0]); xp=(1+u*tau)*e
    return np.array([G*(xp[0]+v*tau),xp[1]])
def weak_thrust(mu,u,v=0.2):
    def emit(s,e):
        X0=boostpos(v,u,e); w=relsum(v,u,e); return X0+w*s, w
    return weak_T(np.zeros(2),np.array([v,0]),emit,0.0,mu,-10.0)
def gr_thrust(mu,u,v=0.2):
    ems=[];hz=0
    for t in th:
        e=np.array([np.cos(t),np.sin(t)]); at,h=thrust(mu,boostpos(v,u,e),relsum(v,u,e),0.0,-10.0); ems.append(at); hz+=h
    Xo=np.zeros(2); Uo=fourvel(mu,*polar(Xo),np.array([v,0]))
    r,bad=observe(mu,Xo,Uo,ems,0.0,-10.0); return r,bad,hz,None
def gr_coast(mu,u,v0,to=3.0):
    ems=[];hz=0;minr=9
    for t in th:
        e=np.array([np.cos(t),np.sin(t)]); at,h=geodesic(mu,e,relsum(v0,u,e),0.0,to); ems.append(at); hz+=h
    oat,h=geodesic(mu,np.zeros(2),np.array([v0,0]),0.0,to); Xo,Uo=oat(to)
    vo=static_vel(mu,*polar(Xo),Uo)
    r,bad=observe(mu,Xo,Uo,ems,to,0.0); return r,bad,hz,(Xo,vo)
def job(c):
    kind,mu,u,v0=c
    try:
        if kind=='thrust': w=weak_thrust(mu,u); g=gr_thrust(mu,u)
        else: w=weak_coast(mu,u,v0); g=gr_coast(mu,u,v0)
    except Exception as ex: return c,None,('error',repr(ex))
    return c,w,g
if __name__=='__main__':
    C=[('thrust',mu,u,0.2) for mu in (0,0.05,0.1,0.2,0.3) for u in (0,0.1)]
    C+=[('coast',mu,u,v0) for v0 in (0.1,0.3,0.5) for u in (0,v0/2,0.9999) for mu in (0,0.02,0.05,0.1)]
    with mp.Pool(8) as P:
        R={}
        for c,w,g in P.imap_unordered(job,C):
            R[c]=(w,g); print(c, None if w is None else "weak %.4f %.4f %.4f"%stats(w)[::3][:1]+(stats(w)[1],stats(w)[3]) if False else '', flush=True)
            pickle.dump(R,open('results.pkl','wb'))
