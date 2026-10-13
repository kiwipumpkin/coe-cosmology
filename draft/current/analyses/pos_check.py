import json, numpy as np
exec(open("scan2.py").read().split("out={")[0])
S=json.load(open("scan.json")); KG=7.7e51/0.0039
# sign-degeneracy test: flip all amplitudes of N=5 sym
p=np.array(S["sym"]["5"]["p"]); sg=S["sym"]["5"]["singles"]; q=p.copy(); q[1:-1:2]*=-1; q[-1]*=-1
print("N5 sym singles",sg,"rms",rms(R-sym_model(p,2,sg)),"flipped-sign rms",rms(R-sym_model(q,2,sg)))
print("corr(model, -model w/ mirror 180):",np.corrcoef(field(0.3,1.7),field(0.3+np.pi,1.7))[0,1])
# variant without monopole removal in field
def field_nm(alpha,D):
    x,y=np.cos(th),np.sin(th); X,Y=D*np.cos(alpha),D*np.sin(alpha); dx,dy=X-x,Y-y
    return (dx*x+dy*y)/(dx*dx+dy*dy)**1.5-(X/D**3*x+Y/D**3*y)
res={}
def run(fam,N,pos=True,fld=field,nrest=4):
    best=None; npair=N//2
    for singles in ([[]] if fam=="gen" or N%2==0 else [[phi],[phi+180]]):
        def mod(p):
            D=p[0]; o=np.zeros_like(th); i=1
            if fam=="sym":
                for _ in range(npair): o+=p[i]*(fld(np.radians(phi+p[i+1]),D)+fld(np.radians(phi-p[i+1]),D)); i+=2
                for s in singles: o+=p[i]*fld(np.radians(s),D); i+=1
            else:
                for j in range(N): o+=p[1+2*j]*fld(np.radians(p[2+2*j]),D)
            return o
        na=npair if fam=="sym" else N
        for D0 in [1.2,1.7,2.5,4]:
            for s in range(nrest):
                g_=np.random.default_rng(s)
                if fam=="sym":
                    p0=[D0]+sum([[abs(g_.normal(0,4e-3)),g_.uniform(5,175)] for _ in range(npair)],[])+[1e-3]*len(singles)
                    lb=[1.05]+[0 if pos else -np.inf,0]*npair+[0 if pos else -np.inf]*len(singles); ub=[20]+[np.inf,180]*npair+[np.inf]*len(singles)
                else:
                    p0=[D0]+sum([[abs(g_.normal(0,4e-3)),g_.uniform(0,360)] for _ in range(N)],[]); lb=[1.05]+[0,-360]*N; ub=[20]+[np.inf,720]*N
                f=least_squares(lambda p:mod(p)-R,p0,bounds=(lb,ub),max_nfev=3000)
                if best is None or f.cost<best[0].cost: best=(f,singles,mod)
    f,singles,mod=best; left=R-mod(f.x)
    amps=(list(f.x[1:1+2*npair:2])+list(f.x[1+2*npair:])) if fam=="sym" else list(f.x[1::2])
    r=dict(rms=rms(left),D_Gpc=float(f.x[0]*14.3),masses_kg=[float(a*KG) for a in amps],modes=modes(left),p=f.x.tolist(),singles=singles)
    print(fam,N,"pos" if pos else "free",f"{r['rms']:.3e}",f"D={r['D_Gpc']:.1f}Gpc",[f"{m:.1e}" for m in r['masses_kg']],r['modes'],flush=True); return r
res["sym5_nomonopole_free"]=run("sym",5,pos=False,fld=field_nm)
res["sym5_nomonopole_pos"]=run("sym",5,pos=True,fld=field_nm)
for N in range(4,9):
    res[f"sym{N}_pos"]=run("sym",N); json.dump(res,open("pos_check.json","w"),indent=1)
for N in range(4,9):
    res[f"gen{N}_pos"]=run("gen",N,nrest=3); json.dump(res,open("pos_check.json","w"),indent=1)
