import sys, os, json, time, numpy as np
OD=os.path.dirname(os.path.abspath(__file__)); SEQ=os.path.dirname(OD); sys.path.insert(0,SEQ)
from seq_core import *
DELTA_C=100.0
def clu_R(M,x): return (3*np.abs(M)/(4*np.pi*DELTA_C*rho_gpc(np.atleast_2d(x))*GPC**3))**(1/3)
def unit(v): return v/np.linalg.norm(v,axis=-1,keepdims=True)
lf=json.load(open(os.path.join(SEQ,"local_fit.json"))); rows=json.load(open(os.path.join(OD,"pairs.json")))
print("same side of ring plane:",np.mean([np.sign(r["void_z"])==np.sign(r["clu_z"]) for r in rows]),
      " cluster closer to plane:",np.mean([abs(r["clu_z"])<abs(r["void_z"]) for r in rows]))

def direction(rule,xv,rng):
    rhat=unit(xv); chat=unit(C_G-xv); phat=unit(np.cross(np.array([0,0,1.0]),xv)); zs=-np.sign(xv[...,2:3]); zs[zs==0]=1
    if rule=="toward_us": return -rhat
    if rule=="away_from_us": return rhat
    if rule=="toward_CoE": return chat
    if rule=="away_from_CoE": return -chat
    if rule=="tangential_+": return phat
    if rule=="tangential_-": return -phat
    if rule=="toward_ring_plane": return np.broadcast_to(np.array([0,0,1.0]),xv.shape)*zs
    if rule=="away_from_ring_plane": return -np.broadcast_to(np.array([0,0,1.0]),xv.shape)*zs
    if rule=="random": return unit(rng.normal(size=xv.shape))
RULES=["toward_us","away_from_us","toward_CoE","away_from_CoE","tangential","toward_ring_plane","away_from_ring_plane","random"]

# ---- (1) keep fitted voids, move each cluster to the rule direction at the SAME offset length (no refit)
def fixed_eval():
    out={}; rng=np.random.default_rng(0)
    for run in lf["runs"]:
        t=run["traj"]; vs=[p for p in t if p["kind"]=="void"]; cs=[p for p in t if p["kind"]=="cluster"]
        XV=np.array([p["x"] for p in vs]); MV=np.array([p["M"] for p in vs]); RVv=np.array([p["R"] for p in vs])
        XC=np.array([p["x"] for p in cs]); L=np.linalg.norm(XC-XV,axis=1)[:,None]
        wv=pattern(XV,MV); res={"as fitted":float(ldiff(wv+pattern(XC,-MV)))}
        for rule in RULES:
            if rule=="random":
                res["random (median of 200)"]=float(np.median([ldiff(wv+pattern(XV+direction("random",XV,rng)*L,-MV)) for _ in range(200)])); continue
            if rule=="tangential":
                res[rule]=float(min(ldiff(wv+pattern(XV+direction("tangential_+",XV,rng)*L,-MV)),ldiff(wv+pattern(XV+direction("tangential_-",XV,rng)*L,-MV)))); continue
            res[rule]=float(ldiff(wv+pattern(XV+direction(rule,XV,rng)*L,-MV)))
        res["voids only (clusters removed)"]=float(ldiff(wv)); res["clusters only (voids removed)"]=float(ldiff(pattern(XC,-MV)))
        out[run["seed"]]=res; print("seed",run["seed"],{k:round(v,1) for k,v in res.items()},flush=True)
    return out

# ---- (2) refit: greedy pairs, cluster forced to sit in the rule direction from its void, offset L in {rim, 3Rv, 5Rv}
def refit(rule,seed,nbank=2500,radii=(1.5,2.25,3.0,3.75),mult=("rim",3.0,5.0),npairs=24):
    rng=np.random.default_rng(seed); pk=Packer(h=2.0); w=np.zeros(360)
    Pv=rand_in_shell(rng,RMIN,LOCAL_RMAX,nbank); Fv=unit_pattern_fast(Pv); rho=rho_gpc(Pv)
    cands=[]
    for Rv in radii:
        okv=region_ok(Pv,Rv,rmax=LOCAL_RMAX); Mv=DELTA_V*rho*vol(Rv)*GPC**3; Rc=clu_R(Mv,Pv)
        dirs=[direction(rule,Pv,rng)] if rule not in ("tangential","random") else ([direction("tangential_+",Pv,rng),direction("tangential_-",Pv,rng)] if rule=="tangential" else [direction("random",Pv,rng)])
        for d in dirs:
            for m in mult:
                L=(Rv+Rc+1.5*GAP) if m=="rim" else np.maximum(m*Rv,Rv+Rc+1.5*GAP)
                Xc=Pv+d*L[:,None] if np.ndim(L) else Pv+d*L
                ok=okv&region_ok(Xc,Rc,rmax=LOCAL_RMAX)
                cands.append((Rv,Mv,Rc,Xc,ok,Mv[:,None]*(Fv-unit_pattern_fast(Xc))))
    for k in range(npairs):
        res=R-w; base=res@res; best=None
        for Rv,Mv,Rc,Xc,ok,F in cands:
            S=np.where(ok,((res[None]-F)**2).sum(1),np.inf)
            for i in np.argsort(S)[:200]:
                if not np.isfinite(S[i]) or S[i]>=base: break
                if pk.free(Pv[i],Rv) and pk.free(Xc[i],Rc[i]) and np.linalg.norm(Pv[i]-Xc[i])>=Rv+Rc[i]+GAP:
                    if best is None or S[i]<best[0]: best=(S[i],Pv[i],Rv,Mv[i],Xc[i],Rc[i],F[i])
                    break
        if best is None: break
        _,xv,Rv,Mv_,xc,Rc_,f=best; pk.add(xv,Rv,Mv_,"void"); pk.add(xc,Rc_,-Mv_,"clu"); w=w+f
    X,r,_,_=pk.arrays()
    return float(ldiff(w)),len(r)//2,float(check_no_overlap(X,r,0.0)) if len(r)>1 else None

if __name__=="__main__":
    if sys.argv[1]=="fixed": json.dump(fixed_eval(),open(os.path.join(OD,"fixed_voids_rules.json"),"w"),indent=1)
    else:
        rule=sys.argv[2]; out=[]
        for s in (0,1):
            L,n,g=refit(rule,s); out.append(dict(seed=s,left=L,pairs=n,min_gap=g)); print(rule,s,round(L,2),n,g,flush=True)
        json.dump(out,open(os.path.join(OD,f"refit_{rule}.json"),"w"))
