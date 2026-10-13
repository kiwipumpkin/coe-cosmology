import json, numpy as np
from scipy.optimize import least_squares
src=open("cmb-cancel.py").read().split('if __name__')[0]; g={}; exec(src,g)
P=json.load(open("cmb-cancel-story.json"))["params"]; phi=P[4]
ELL=g["ELL"]; R=g["resid"](g["G_FULL"],g["S_model"](*P)[0]); R=R-R.mean()
th=np.radians(ELL); rms=lambda x: float(np.sqrt(np.mean(x**2)))
def field(alpha,D):  # identical to clump-fit.py
    x,y=np.cos(th),np.sin(th); X,Y=D*np.cos(alpha),D*np.sin(alpha)
    dx,dy=X-x,Y-y; r3=(dx*dx+dy*dy)**1.5
    ar=(dx*x+dy*y)/r3 - (X/D**3*x+Y/D**3*y); return ar-ar.mean()
def sym_model(p,npair,singles):
    D=p[0]; out=np.zeros_like(th); i=1
    for _ in range(npair):
        k,off=p[i],p[i+1]; i+=2
        out+=k*(field(np.radians(phi+off),D)+field(np.radians(phi-off),D))
    for s in singles: out+=p[i]*field(np.radians(s),D); i+=1
    return out
def gen_model(p,N):
    D=p[0]; out=np.zeros_like(th)
    for j in range(N): out+=p[1+2*j]*field(np.radians(p[2+2*j]),D)
    return out
def modes(left):
    v=np.mean(left**2); r=[]
    for m in range(1,16):
        A=np.hypot(2*np.mean(left*np.cos(m*th)),2*np.mean(left*np.sin(m*th))); r.append((m,0.5*A*A/v))
    r.sort(key=lambda t:-t[1]); return [f"m={m}({100*f:.0f}%)" for m,f in r[:3]]
def fit(fun,p0s,lb,ub):
    best=None
    for p0 in p0s:
        try: f=least_squares(lambda p:fun(p)-R,p0,bounds=(lb,ub),max_nfev=4000)
        except Exception: continue
        if best is None or f.cost<best.cost: best=f
    return best
out={"n_samples":len(th),"rms_post_SR":rms(R),"sym":{},"gen":{}}
print("samples",len(th),"post-SR rms",rms(R))
for N in range(4,9):
    # A: original symmetric family (pairs about the axis, shared D; odd N adds on-axis clump)
    npair=N//2; best=None
    for singles in ([[]] if N%2==0 else [[phi],[phi+180]]):
        p0s=[]
        for D0 in [1.2,1.5,2,3,4,8]:
            for s in range(8):
                rng=np.random.default_rng(s)
                p0s.append([D0]+sum([[rng.normal(0,1e-2),rng.uniform(5,175)] for _ in range(npair)],[])+[0.0]*len(singles))
        lb=[1.05]+[-np.inf,0]*npair+[-np.inf]*len(singles); ub=[20]+[np.inf,180]*npair+[np.inf]*len(singles)
        f=fit(lambda p:sym_model(p,npair,singles),p0s,lb,ub)
        if best is None or f.cost<best[0].cost: best=(f,singles)
    f,singles=best; left=R-sym_model(f.x,npair,singles)
    ks=np.abs(np.r_[f.x[1:1+2*npair:2],f.x[1+2*npair:]])
    out["sym"][N]=dict(singles=[float(s) for s in singles],rms=rms(left),nparam=len(f.x),D=float(f.x[0]),max_abs_amp=float(ks.max()),modes=modes(left),p=f.x.tolist())
    # B: general family (free angle + amplitude per clump, shared D)
    p0s=[]
    for D0 in [1.2,1.5,2,3,5]:
        for s in range(10):
            rng=np.random.default_rng(100+s)
            p0s.append([D0]+sum([[rng.normal(0,1e-2),rng.uniform(0,360)] for _ in range(N)],[]))
    f=fit(lambda p:gen_model(p,N),p0s,[1.05]+[-np.inf,-360]*N,[20]+[np.inf,720]*N)
    left=R-gen_model(f.x,N)
    out["gen"][N]=dict(rms=rms(left),nparam=len(f.x),D=float(f.x[0]),max_abs_amp=float(np.abs(f.x[1::2]).max()),modes=modes(left),p=f.x.tolist())
    for k in ("sym","gen"):
        r=out[k][N]; json.dump(out,open("scan.json","w"),indent=1); print(N,k,f"{r['rms']:.3e}",r['nparam'],f"D={r['D']:.2f}",f"maxk={r['max_abs_amp']:.2e}",r['modes'],flush=True)
json.dump(out,open("scan.json","w"),indent=1)
