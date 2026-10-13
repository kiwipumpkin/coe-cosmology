"""Scenario A. (1) Fit a local arrangement (14.3-30 Gpc) by adding voids and clusters ONE AT A TIME (void, then its
matching cluster), every sphere non-overlapping, each pair chosen to lower the leftover the most. (2) Freeze it.
(3) Fill the outer part of the surround region (30-42.9 Gpc from us, CoE interior excluded) with compensated
void+rim-cluster units placed one at a time, in 4 ways: empty / random / foam-like (grown touching neighbours) / lattice.
usage: python3 scenario_a.py local          -> local_fit.json, local_frozen.npz
       python3 scenario_a.py outer <cfg> <Rv> <nreal> <seed>   -> outer_<cfg>_Rv<Rv>.json  (cfg: random|foam|lattice)"""
import sys, os, json, time, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from seq_core import *
DELTA_C=100.0                       # cluster contrast: keeps clusters ~0.5% of volume vs voids' share (PACK)
def clu_R(M,x): return (3*abs(M)/(4*np.pi*DELTA_C*rho_gpc(np.asarray(x))*GPC**3))**(1/3)

# ---------------------------------------------------------------- local fit
def local_fit(seed=0,nbank=6000,radii=(1.5,2.25,3.0,3.75),max_pairs=24,top=1500):
    rng=np.random.default_rng(seed); pk=Packer(h=2.0); w=np.zeros(360)
    # candidate bank (positions in near region), fixed for this seed
    Pv=rand_in_shell(rng,RMIN,LOCAL_RMAX,nbank); Pc=rand_in_shell(rng,RMIN,LOCAL_RMAX,nbank)
    Fv=unit_pattern_fast(Pv); Fc=unit_pattern_fast(Pc)
    rhoV=rho_gpc(Pv)
    traj=[dict(step=0,kind="start",left=float(ldiff(w)),amps=(amps(R)/FLOOR)[:6].tolist())]
    for k in range(max_pairs):
        res=R-w; base=float(res@res); best=None
        for Rv in radii:
            okv=region_ok(Pv,Rv,rmax=LOCAL_RMAX); 
            if not okv.any(): continue
            Mv=DELTA_V*rhoV*vol(Rv)*GPC**3                 # (nbank,) negative
            # pair effect: Mv_i*Fv_i + (-Mv_i)*Fc_j ; choose top voids & clusters by linearised gain, then exact pair search
            gv=np.where(okv,-(Mv*(Fv@res)),-np.inf)         # want Mv*Fv aligned with res -> gain = Mv*(Fv.res)
            iv=np.argsort(-gv)[:top]; iv=iv[np.isfinite(gv[iv])]
            for sgn_sel in (1,):
                pass
            A=Mv[iv,None]*Fv[iv]                             # void patterns (top,360)
            # for cluster we need mass -Mv_i; pattern = -Mv_i*Fc_j
            # |res - A_i + Mv_i Fc_j|^2 = |res-A_i|^2 + 2 Mv_i (res-A_i).Fc_j + Mv_i^2 |Fc_j|^2
            RA=res[None]-A
            t1=(RA*RA).sum(1)                                # (top,)
            t2=2*Mv[iv,None]*(RA@Fc.T)                       # (top,nbank)
            t3=(Mv[iv,None]**2)*((Fc*Fc).sum(1))[None]
            S=t1[:,None]+t2+t3
            order=np.argsort(S,axis=None)[:4000]
            for o in order:
                i,j=divmod(int(o),S.shape[1]); vi=iv[i]
                if S[i,j]>=base: break
                xv=Pv[vi]; xc=Pc[j]; Mvi=Mv[vi]; Rc=clu_R(Mvi,xc)
                if not pk.free(xv,Rv): continue
                if not region_ok(xc,Rc,rmax=LOCAL_RMAX)[0]: continue
                if not pk.free(xc,Rc): continue
                if np.linalg.norm(xv-xc)<Rv+Rc+GAP: continue
                if best is None or S[i,j]<best[0]: best=(S[i,j],xv,Rv,Mvi,xc,Rc)
                break
        if best is None or best[0]>=base*(1-1e-4): break
        _,xv,Rv,Mvi,xc,Rc=best
        pk.add(xv,Rv,Mvi,"void"); w=w+pattern(xv[None],np.array([Mvi]))
        traj.append(dict(step=len(traj),kind="void",left=float(ldiff(w)),R=Rv,x=xv.tolist(),M=Mvi))
        pk.add(xc,Rc,-Mvi,"clu"); w=w+pattern(xc[None],np.array([-Mvi]))
        traj.append(dict(step=len(traj),kind="cluster",left=float(ldiff(w)),R=float(Rc),x=xc.tolist(),M=-Mvi,amps=(amps(R-w)/FLOOR)[:6].tolist()))
        print(f"seed {seed} pair {k+1}: void R={Rv} Gpc  left={traj[-1]['left']:.1f}x",flush=True)
    X,r,M,kind=pk.arrays()
    return dict(seed=seed,X=X,r=r,M=M,kind=kind,w=w,traj=traj,left=float(ldiff(w)),margin=check_no_overlap(X,r,0.0) if len(X) else np.inf)

# ---------------------------------------------------------------- outer fills
OUT_R0=LOCAL_RMAX
def unit_spheres(x,Rv,rng,nrim=6):
    """void at x + nrim clusters on its rim (random orientation); each cluster carries 1/nrim of the deficit."""
    Mv=float(void_mass(x[None],Rv)[0]); Mc=-Mv/nrim; Rc=clu_R(Mc,x)
    Q,_=np.linalg.qr(rng.normal(size=(3,3)))
    dirs=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)@Q.T
    xc=x+dirs*(Rv+Rc+GAP*1.5)
    return [(x,Rv,Mv,"void")]+[(c,Rc,Mc,"clu") for c in xc]

RIM="random"
def _free_planned(c,r,plan):
    for p,pr,_,_ in plan:
        d=p-c
        if d@d<(r+pr+GAP)**2: return False
    return True
def try_unit(pk,x,Rv,rng,nrim=6):
    """void first, then its clusters ONE AT A TIME at free spots on its rim (random directions -> irregular junctions);
    RIM='octa' gives the perfectly symmetric rim (invisible-unit idealisation)."""
    if not region_ok(x,Rv,rmin=OUT_R0)[0] or not pk.free(x,Rv): return False
    if RIM=="octa":
        sp=unit_spheres(x,Rv,rng,nrim)
        for c,r,_,_ in sp:
            if not region_ok(c,r,rmin=OUT_R0)[0] or not pk.free(c,r): return False
        for c,r,M,k in sp: pk.add(c,r,M,k)
        return True
    Mv=float(void_mass(x[None],Rv)[0]); Mc=-Mv/nrim; Rc=clu_R(Mc,x); plan=[(x,Rv,Mv,"void")]
    for _ in range(nrim):
        ok=False
        for _t in range(25):
            c=x+rand_dir(rng)*(Rv+Rc+GAP*1.5)
            if region_ok(c,Rc,rmin=OUT_R0)[0] and pk.free(c,Rc) and _free_planned(c,Rc,plan): ok=True; break
        if not ok: return False
        plan.append((c,Rc,Mc,"clu"))
    for c,r,M,k in plan: pk.add(c,r,M,k)
    return True

def fill_random(pk,Rv,rng,max_fail=4000):
    n=0; fail=0
    while fail<max_fail:
        x=rand_in_shell(rng,OUT_R0,S_MAX_G,1)[0]
        if try_unit(pk,x,Rv,rng): n+=1; fail=0
        else: fail+=1
    return n

def fill_foam(pk,Rv,rng,max_fail=4000):
    """grow outward: each new unit sits right next to an already-placed void (shared rim), like a foam of cells."""
    Rc=clu_R(float(void_mass(np.array([[0,35,0.]]),Rv)[0])/6,np.array([[0,35,0.]]))[0]
    step=2*Rv+2*Rc+3*GAP                              # touching cell spacing (rims shared)
    fronts=[]; n=0; fail=0
    # seed: 8 random units
    while len(fronts)<8 and fail<max_fail:
        x=rand_in_shell(rng,OUT_R0,S_MAX_G,1)[0]
        if try_unit(pk,x,Rv,rng): fronts.append(x); n+=1
        else: fail+=1
    fail=0
    while fronts:
        j=rng.integers(len(fronts)); p=fronts[j]; ok=False
        for _ in range(24):
            x=p+rand_dir(rng)*step*rng.uniform(1.0,1.15)
            if try_unit(pk,x,Rv,rng): fronts.append(x); n+=1; ok=True; break
        if not ok: fronts.pop(j)
    return n

def fill_lattice(pk,Rv,rng):
    """simple cubic lattice of units, spacing just large enough for no overlap; random rotation and offset."""
    Rc=clu_R(float(void_mass(np.array([[0,35,0.]]),Rv)[0])/6,np.array([[0,35,0.]]))[0]
    a=2*(Rv+2*Rc+2*GAP)*1.02
    Q,_=np.linalg.qr(rng.normal(size=(3,3))); off=rng.uniform(0,a,3)
    n=int(np.ceil(S_MAX_G/a))+1; g=np.arange(-n,n+1)
    G=np.stack(np.meshgrid(g,g,g,indexing="ij"),-1).reshape(-1,3)*a+off
    G=G@Q.T; s=np.linalg.norm(G,axis=1); G=G[(s>OUT_R0-a)&(s<S_MAX_G+a)]
    G=G[np.argsort(np.linalg.norm(G,axis=1))]
    cnt=0
    for x in G:
        # lattice orientation for rims too (aligned with lattice axes) -> use fixed Q
        Mv=float(void_mass(x[None],Rv)[0]); Mc=-Mv/6; Rc_=clu_R(Mc,x)
        dirs=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)@Q.T
        # rotate rims 45deg about one axis so neighbours' rims don't collide head-on
        sp=[(x,Rv,Mv,"void")]+[(x+d*(Rv+Rc_+GAP*1.5),Rc_,Mc,"clu") for d in dirs]
        if all(region_ok(c,r,rmin=OUT_R0)[0] and pk.free(c,r) for c,r,_,_ in sp):
            for c,r,M,k in sp: pk.add(c,r,M,k)
            cnt+=1
    return cnt

def fill_pairs(pk,Rv,rng,max_fail=4000,radii=(1.5,2.25,3.0,3.75),n_pairs=48):
    """control with the SAME kind of objects as the local fix: big voids and single heavy clusters (mass = that void's deficit),
    placed one at a time at random spots (cluster not tied to its void). n_pairs per layout (24 = as many as the local fix)."""
    n=0; fail=0
    while n<n_pairs and fail<max_fail:
        Rv_=radii[rng.integers(len(radii))]; x=rand_in_shell(rng,OUT_R0,S_MAX_G,1)[0]
        if not (region_ok(x,Rv_,rmin=OUT_R0)[0] and pk.free(x,Rv_)): fail+=1; continue
        Mv=float(void_mass(x[None],Rv_)[0])
        for _ in range(200):
            c=rand_in_shell(rng,OUT_R0,S_MAX_G,1)[0]; Rc=clu_R(Mv,c)
            if region_ok(c,Rc,rmin=OUT_R0)[0] and pk.free(c,Rc) and np.linalg.norm(c-x)>=Rv_+Rc+GAP: break
        else: fail+=1; continue
        pk.add(x,Rv_,Mv,"void"); pk.add(c,Rc,-Mv,"clu"); n+=1
    return n
def fill_lattice_octa(pk,Rv,rng):
    return fill_lattice(pk,Rv,rng)
FILLS=dict(random=fill_random,foam=fill_foam,lattice=fill_lattice,pairs=fill_pairs)

def load_local():
    z=np.load(os.path.join(HERE,"local_frozen.npz"),allow_pickle=True); return z

def outer(cfg,Rv,nreal,seed):
    z=load_local(); w_loc=z["w"]; L0=float(ldiff(w_loc))
    fn=os.path.join(HERE,f"outer_{cfg}_Rv{Rv:g}_s{seed}.json")
    rows=json.load(open(fn))["rows"] if os.path.exists(fn) else []
    t0=time.time()
    for i in range(len(rows),nreal):
        rng=np.random.default_rng([seed,i]); pk=Packer(h=max(2.0,2*Rv))
        for x,r,M,k in zip(z["X"],z["r"],z["M"],z["kind"]): pk.add(x,r,M,str(k))
        n0=len(pk.X)
        nu=FILLS[cfg](pk,Rv,rng)
        X,r,M,kind=pk.arrays(); Xo,Mo=X[n0:],M[n0:]
        w_out=pattern(Xo,Mo) if len(Mo) else np.zeros(360)
        wR0=ldiff(w_out)                                  # same outer fill with NO local fix (vs 237)
        L=ldiff(w_loc+w_out)
        void_frac=float(np.sum(vol(r[n0:][kind[n0:]=="void"]))/180678.6)
        margin=check_no_overlap(X,r,0.0)
        rows.append(dict(i=i,units=int(nu),left=float(L),delta=float(L-L0),left_without_local=float(wR0),
                         outer_rms=float(rms(w_out)/FLOOR),void_vol_frac=void_frac,net_kg=float(Mo.sum()),min_gap=float(margin),
                         amps=(amps(R-w_loc-w_out)/FLOOR)[:6].tolist()))
        print(f"{cfg} Rv={Rv} {i+1}/{nreal}: units={nu} voidfrac={void_frac:.3f} left={L:.2f} (frozen {L0:.2f}) outer_rms={rms(w_out)/FLOOR:.2f} gap={margin:.3f} t={time.time()-t0:.0f}s",flush=True)
        json.dump(dict(cfg=cfg,Rv=Rv,frozen_left=L0,rows=rows),open(fn,"w"),indent=1)

if __name__=="__main__":
    if sys.argv[1]=="local":
        runs=[]
        for sd in range(int(sys.argv[2]) if len(sys.argv)>2 else 4):
            r=local_fit(seed=sd); runs.append(r); print(f"seed {sd}: final {r['left']:.2f}x  objects={len(r['M'])} min_gap={r['margin']:.3f}",flush=True)
        b=min(runs,key=lambda r:r["left"])
        np.savez(os.path.join(HERE,"local_frozen.npz"),X=b["X"],r=b["r"],M=b["M"],kind=b["kind"],w=b["w"])
        json.dump(dict(best_seed=b["seed"],runs=[dict(seed=r["seed"],left=r["left"],n_obj=len(r["M"]),min_gap=r["margin"],traj=r["traj"],
                       final_amps=(amps(R-r["w"])/FLOOR)[:6].tolist()) for r in runs],
                       start_amps=(amps(R)/FLOOR)[:6].tolist(),four_clump_amps=(amps(R4)/FLOOR)[:6].tolist()),
                  open(os.path.join(HERE,"local_fit.json"),"w"),indent=1)
    else:
        global_rim=sys.argv[6] if len(sys.argv)>6 else "random"
        RIM=global_rim
        outer(sys.argv[2],float(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]))
