import sys,pickle; sys.path.insert(0,'.'); from boost import *
mu=0.1
def t180(v0):
    g=1/np.sqrt(1-v0*v0); ots,oX=g*v0,np.array([g*v0*v0,0])
    def ok(to):
        oat,_=geod2(mu,oX,np.array([v0,0]),ots,to); Xo,_=oat(to); so,L=ray(mu,Xo,0.0,to,0.0)
        return cart(*so.sol(0.0)[:2])[0] < -1/g
    lo,hi=max(ots,2*g*v0)+0.01,max(3.0,3/(1-v0))
    while not ok(hi): hi*=1.5
    for _ in range(30):
        m=(lo+hi)/2
        if ok(m): hi=m
        else: lo=m
    return hi
def job(v0):
    t=t180(v0); step=0.1 if t<5 else (0.5 if t<30 else 1.0); to=np.ceil(t/step-1e-9)*step
    while True:
        res={}
        for u in (0,0.9999):
            for m in (0.0,mu):
                g=gr_boost(m,u,v0,to)[0]; res[(m,u)]=g
        if not any(np.isnan(r).any() for r in res.values()): return v0,to,t,res
        to+=step
if __name__=='__main__':
    V=(0.01,0.05,0.1,0.2,0.3,0.5,0.7,0.9,0.95,0.99)
    with mp.Pool(8) as P:
        out={}
        for v0,to,t,res in P.imap_unordered(job,V):
            out[v0]=(to,res); pickle.dump(out,open('sweep_static.pkl','wb'))
            s=lambda r:"%.4f %.4f %.4f"%(stats(norm(r))[0],stats(norm(r))[1],stats(norm(r))[3])
            print(v0,"view t",to,"(180 deg limit %.3f)"%t," | ".join(f"m{m} u{u}: "+s(r) for (m,u),r in res.items()),flush=True)
