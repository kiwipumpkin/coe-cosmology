import json, numpy as np
exec(open("scan2.py").read().split("out={")[0])
S=json.load(open("scan.json")); KG=7.7e51/0.0039; out={}
def sym_seed(N):
    p=S["sym"][str(N)]["p"]; sg=S["sym"][str(N)]["singles"]; s=[p[0]]; i=1
    for _ in range(N//2): s+=[abs(p[i]),phi+p[i+1],abs(p[i]),phi-p[i+1]]; i+=2
    for a in sg: s+=[abs(p[i]),a]; i+=1
    return s
for Dmax,tag in [(20,"Dle20"),(4,"Dle4")]:
  prev=None
  for N in range(4,9):
    seeds=[sym_seed(N)]+[[min(sym_seed(N)[0],Dmax-0.01)]+sym_seed(N)[1:]]
    if prev is not None:
        for a in np.arange(0,360,30): seeds.append(list(prev)+[1e-3,a])
    for s in range(24):
        g_=np.random.default_rng(500+s); seeds.append([g_.uniform(1.2,min(Dmax,5))]+sum([[g_.uniform(1e-3,8e-3),g_.uniform(0,360)] for _ in range(N)],[]))
    best=None
    for p0 in seeds:
        p0=np.array(p0,float); p0[0]=np.clip(p0[0],1.06,Dmax-1e-3); p0[1::2]=np.maximum(p0[1::2],1e-6)
        f=least_squares(lambda p:gen_model(p,N)-R,p0,bounds=([1.05]+[0,-720]*N,[Dmax]+[np.inf,1080]*N),max_nfev=3000)
        if best is None or f.cost<best.cost: best=f
    p=best.x; prev=p; left=R-gen_model(p,N)
    r=dict(rms=rms(left),D=float(p[0]),D_Gpc=float(p[0]*14.3),angles_deg=[float(a%360) for a in p[2::2]],masses_kg=[float(k*KG) for k in p[1::2]],modes=modes(left),p=p.tolist())
    out[f"{tag}_N{N}"]=r; json.dump(out,open("pos_free.json","w"),indent=1)
    print(tag,N,f"{r['rms']:.3e}",f"D={r['D_Gpc']:.1f}",[f"{a:.0f}" for a in r['angles_deg']],[f"{m:.1e}" for m in r['masses_kg']],r['modes'],flush=True)
