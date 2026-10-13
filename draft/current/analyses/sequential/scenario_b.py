"""Scenario B: grow structure from empty, ONE object at a time (a void or a cluster), each placement conditioned on what is
already there, strict non-overlap. Region: 14.3-42.9 Gpc from us, CoE interior excluded (same surround region as before).
Voids R=Rv (delta=-0.70); clusters delta=+100 (compact, ~0.5% of volume), mass balance kept by clusters paying back voids.
Leftover is recorded after every object.
usage: python3 scenario_b.py <rule> <seed> [Rv]
rules: foam_rim    - new void grows next to an existing void; then its 6 clusters one by one at rim spots touching other voids (junctions)
       alternate   - void, cluster, void, cluster...: one heavy cluster per void (all of its deficit), at a junction across a NEIGHBOURING void
       coe_bias    - like foam_rim but voids grow preferentially on the side away from the CoE, clusters on CoE-facing rims
       random      - same objects/mix as foam_rim, each one dropped at a random free spot (no conditioning) -> unstructured control
       cmb_aware   - like foam_rim but each placement is the best of ~30 allowed spots for the leftover (structure + aim)"""
import sys, os, json, time, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from seq_core import *
DELTA_C=100.0; NRIM=6
def clu_R(M,x): return float((3*abs(M)/(4*np.pi*DELTA_C*rho_gpc(np.asarray(x).reshape(1,3))[0]*GPC**3))**(1/3))
AXC=C_G/np.linalg.norm(C_G)
def ok(pk,x,r): return bool(region_ok(x,r)[0]) and pk.free(x,r)

class Grow:
    def __init__(s,rule,seed,Rv=1.5,max_fail=3000):
        s.rule=rule; s.rng=np.random.default_rng(seed); s.Rv=Rv; s.pk=Packer(h=max(2.0,2*Rv)); s.voids=[]; s.fronts=[]
        s.max_fail=max_fail; s.w=np.zeros(360); s.order=[]   # (x,M) in placement order
        Mv0=float(void_mass(np.array([[0,30,0.]]),Rv)[0]); s.Rc6=clu_R(Mv0/NRIM,[0,30,0]); s.Rc1=clu_R(Mv0,[0,30,0])
    def put(s,x,r,M,k):
        s.pk.add(x,r,M,k); s.order.append((np.array(x),M,k))
        if s.rule=="cmb_aware": s.w=s.w+pattern(np.array(x)[None],np.array([M]))
    def pick(s,cands,Ms):
        """cmb_aware: choose candidate that lowers leftover most; others: first candidate."""
        if s.rule!="cmb_aware" or len(cands)==1: return 0
        F=unit_pattern_fast(np.array(cands))*np.array(Ms)[:,None]; res=R-s.w
        return int(np.argmin(((res[None]-F)**2).sum(1)))
    def void_spot(s,n=30):
        rng=s.rng; out=[]
        for _ in range(400):
            if s.rule=="random" or not s.fronts:
                x=rand_in_shell(rng,RMIN,S_MAX_G,1)[0]
            else:
                j=rng.integers(len(s.fronts)); p,Rp=s.fronts[j]
                x=p+rand_dir(rng)*(Rp+s.Rv+2*s.Rc6+4*GAP)*rng.uniform(1.0,1.12)
            if s.rule=="coe_bias":
                cos=float(np.dot(x/np.linalg.norm(x),-AXC))
                if rng.uniform()>(1+0.8*cos)/1.8: continue
            if ok(s.pk,x,s.Rv): out.append(x)
            if len(out)>=(n if s.rule=="cmb_aware" else 1): break
        return out
    def rim_spots(s,v,Rv,Rc,n=25):
        rng=s.rng; out=[]; sc=[]
        V=np.array([q for q,_ in s.voids]) if s.voids else np.zeros((0,3))
        for _ in range(n*3):
            d=rand_dir(rng)
            if s.rule=="coe_bias" and np.dot(d,AXC)<0: d=-d
            c=v+d*(Rv+Rc+1.5*GAP)
            if s.rule=="random": c=rand_in_shell(rng,RMIN,S_MAX_G,1)[0]
            if ok(s.pk,c,Rc):
                # junction score: how many void rims this cluster touches
                touch=int(np.sum(np.abs(np.linalg.norm(V-c,axis=1)-Rv-Rc)<0.6)) if len(V) else 0
                out.append(c); sc.append(touch)
            if len(out)>=n: break
        return out,sc
    def step_void(s):
        c=s.void_spot()
        if not c: return None
        Mv=[float(void_mass(x[None],s.Rv)[0]) for x in c]; i=s.pick(c,Mv); x=c[i]
        s.put(x,s.Rv,Mv[i],"void"); s.voids.append((x,s.Rv)); s.fronts.append((x,s.Rv)); return x,Mv[i]
    def step_cluster(s,v,M,Rc,around=None):
        center=v if around is None else around
        c,sc=s.rim_spots(center,s.Rv,Rc)
        if not c: return False
        if s.rule=="cmb_aware": i=s.pick(c,[M]*len(c))
        else:
            sc=np.array(sc); best=np.where(sc==sc.max())[0]; i=int(best[0])
        s.put(c[i],Rc,M,"clu"); return True
    def run(s):
        fail=0; t0=time.time()
        while fail<20 and len(s.voids)<getattr(s,'max_voids',10**9):
            r=s.step_void()
            if r is None:
                # front exhausted -> drop random fronts and retry; stop when region jams
                if s.fronts: s.fronts=s.fronts[len(s.fronts)//2:]
                fail+=1; continue
            x,Mv=r; fail=0
            if s.rule=="alternate":
                # one heavy cluster: on the far side of a neighbouring void (junction across a neighbour), else own rim
                nb=[q for q,_ in s.voids[:-1] if np.linalg.norm(q-x)<2*s.Rv+2*s.Rc6+1.0]
                Rc=clu_R(-Mv,x); placed=False
                if nb:
                    q=nb[s.rng.integers(len(nb))]; placed=s.step_cluster(x,-Mv,Rc,around=q)
                if not placed: placed=s.step_cluster(x,-Mv,Rc)
                if not placed: s.undo_void(); continue
            else:
                Rc=clu_R(-Mv/NRIM,x); n=0
                for _ in range(NRIM):
                    if s.step_cluster(x,-Mv/NRIM,Rc): n+=1
                if n<NRIM:   # could not pay back fully -> remove the clusters & void just placed (keeps net zero exactly)
                    s.undo_unit(1+n); continue
            if len(s.voids)%200==0: print(f"{s.rule} voids={len(s.voids)} objs={len(s.order)} t={time.time()-t0:.0f}s",flush=True)
        return s
    def undo_void(s): s.undo_unit(1)
    def undo_unit(s,n):
        # rebuild packer without the last n objects (rare)
        keep=s.order[:-n]; X,r,M,kind=s.pk.arrays(); s.pk=Packer(h=s.pk.h)
        for x,rr,m,k in zip(X[:-n],r[:-n],M[:-n],kind[:-n]): s.pk.add(x,rr,m,str(k))
        s.order=keep; s.voids=s.voids[:-1]; s.fronts=[f for f in s.fronts if not np.array_equal(f[0],X[-n])]
        if s.rule=="cmb_aware": s.w=pattern(np.array([o[0] for o in keep]),np.array([o[1] for o in keep])) if keep else np.zeros(360)

def trajectory(order):
    X=np.array([o[0] for o in order]); M=np.array([o[1] for o in order]); k=np.array([o[2] for o in order])
    w=np.zeros(360); left=np.empty(len(M)); left4=np.empty(len(M))
    for i in range(0,len(M),2000):
        F=unit_pattern_fast(X[i:i+2000])*M[i:i+2000,None]; C=np.cumsum(F,0)+w
        left[i:i+2000]=np.sqrt(np.mean((R[None]-C)**2,1))/FLOOR
        w=C[-1]
    return X,M,k,left,w

if __name__=="__main__":
    rule=sys.argv[1]; seed=int(sys.argv[2]); Rv=float(sys.argv[3]) if len(sys.argv)>3 else 1.5
    t0=time.time(); g=Grow(rule,seed,Rv).run()
    X,M,k,left,w=trajectory(g.order); _,r,_,_=g.pk.arrays()
    gap=check_no_overlap(g.pk.arrays()[0],r,0.0)
    nv=np.cumsum(k=="void")
    # store trajectory: every object for first 2000, then every 10th
    idx=np.r_[np.arange(min(2000,len(M))),np.arange(2000,len(M),10)].astype(int) if len(M)>2000 else np.arange(len(M))
    vf=float(np.sum(k=="void")*vol(Rv)/268960.0)
    out=dict(rule=rule,seed=seed,Rv=Rv,n_obj=int(len(M)),n_void=int((k=="void").sum()),n_clu=int((k=="clu").sum()),
             void_vol_frac=vf,net_kg=float(M.sum()),min_gap_gpc=float(gap),final_left=float(left[-1]),min_left=float(left.min()),
             web_rms=float(rms(w)/FLOOR),final_amps=(amps(R-w)/FLOOR)[:6].tolist(),
             traj_obj=(idx+1).tolist(),traj_void=nv[idx].tolist(),traj_left=left[idx].round(3).tolist(),
             traj_kind=[str(x) for x in k[idx[:400]]],elapsed=time.time()-t0)
    json.dump(out,open(os.path.join(HERE,f"b_{rule}_Rv{Rv:g}_s{seed}.json"),"w"))
    print(f"DONE {rule} seed={seed} voids={out['n_void']} clusters={out['n_clu']} voidfrac={vf:.3f} final={left[-1]:.2f}x min={left.min():.2f}x web_rms={out['web_rms']:.2f} gap={gap:.3f} net={M.sum():.2e} t={out['elapsed']:.0f}s",flush=True)
