"""Global search over orientation (3 Euler angles rel. to the CoE axis) and phase (3) for each lattice type and spacing,
amplitude fitted in closed form (bounded); target = R (237x start) or R4 (four compensated clumps, 88.8x).
run: python3 run_search.py  -> search_results.json"""
import json, sys, time, numpy as np
from multiprocessing import Pool
from scipy.stats import qmc
from scipy.optimize import minimize
exec(open("lattice_core.py").read())
SPACINGS=[5,6,7,8,9,10,12,14,16,18,20,22,25,28.6,32,36,40]
KINDS=["SC_rocksalt","BCC_CsCl","FCC_L10","HCP_AB","shell_rings"]
TARGETS={"alone":R,"with4":R4}

BUFFER=True   # main rule: whole object (delta=+1 cluster / delta=-1 void sphere holding |M|) outside the 14.3 Gpc vision sphere
def unpack(kind,x):
    """x in unit cube -> uniform rotation (ZYZ Euler, beta=arccos(1-2u)) and phase; periodic/folded so Nelder-Mead can roam."""
    u=np.asarray(x,float)
    if kind=="shell_rings": return (0.0,0.0,0.0),tuple(np.mod(u[:3],1.0))
    v=np.mod(u[1],2.0); v=v if v<=1 else 2-v
    return (2*np.pi*np.mod(u[0],1),float(np.arccos(1-2*v)),2*np.pi*np.mod(u[2],1)),tuple(np.mod(u[3:6],1.0))

def evaluate(kind,a,x,target,buffer=None):
    buffer=BUFFER if buffer is None else buffer
    e,ph=unpack(kind,x); P,s=lattice_sites(kind,a,e,ph)
    for it in range(10):      # buffer: drop objects that would poke into the vision sphere at the fitted strength, refit, repeat
        if len(P)==0 or (s>0).sum()==0 or (s<0).sum()==0: return rms(target),None
        w,vs,Vs,M=pattern(P,s,a,kind); A,r=best_amp(target,w,vs)
        if not buffer: break
        keep=np.linalg.norm(P,axis=1)>1.0+radius_gpc(np.abs(M*A),1.0)/RV
        if keep.all(): break
        P,s=P[keep],s[keep]
    else:
        return rms(target),None   # did not settle (never seen in practice)
    return r,dict(P=P,s=s,w=w,vs=vs,Vs=Vs,M=M*A,A=A,euler=e,phase=ph)

def task(args):
    kind,aG,tname=args; a=aG/RV; target=TARGETS[tname]; base=rms(target); t0=time.time()
    dim=3 if kind=="shell_rings" else 6
    nS=512 if kind=="shell_rings" else (1024 if aG<10 else 2048)
    X=qmc.Sobol(dim,seed=int(aG*10)+7).random(nS)
    f=np.array([evaluate(kind,a,x,target)[0] for x in X])
    order=np.argsort(f); best=(f[order[0]],X[order[0]])
    npol=3 if aG<8 else 6; mfev=150 if aG<8 else 400
    for i in order[:npol]:
        r=minimize(lambda x:evaluate(kind,a,x,target)[0],X[i],method="Nelder-Mead",
                   options=dict(maxfev=mfev,xatol=1e-4,fatol=1e-9,initial_simplex=X[i]+0.03*np.vstack([np.zeros(dim),np.eye(dim)])))
        if r.fun<best[0]: best=(r.fun,r.x)
    r,info=evaluate(kind,a,best[1],target)
    Mabs=np.abs(info["M"]); sd=np.linalg.norm(info["P"],axis=1); rc=np.linalg.norm(info["P"]-C,axis=1)
    # buffer check: object (cluster delta=+1 / void delta=-1 sphere holding |M|) entirely outside the 14.3 Gpc vision sphere?
    robj=radius_gpc(Mabs,1.0)/RV; inside=int((sd<1+robj).sum())
    out=dict(rule="buffer" if BUFFER else "bare",kind=kind,spacing_gpc=aG,target=tname,base=base/FLOOR,best=r/FLOOR,sobol_median=float(np.median(f))/FLOOR,sobol_p5=float(np.percentile(f,5))/FLOOR,
             n_sobol=nS,x=list(map(float,best[1])),euler_deg=[float(np.degrees(v)) for v in info["euler"]],phase=list(map(float,info["phase"])),
             amp=info["A"],amp_max=amp_max(info["vs"]),void_rescale=info["vs"],n_sites=int(len(info["P"])),n_clusters=int((info["s"]>0).sum()),n_voids=int((info["s"]<0).sum()),
             cell_mass_kg=float(RHO*info["Vs"]),M_obj_median_kg=float(np.median(Mabs)),M_obj_max_kg=float(Mabs.max()),
             net_mass_kg=float(info["M"].sum()),min_dist_from_us_gpc=float(sd.min()*RV),n_within_buffer=inside,
             frac_beyond_Dpd=float(np.mean(rc>B_SH)),amps_over_floor=(amps(target-info["A"]*info["w"])/FLOOR).tolist(),secs=time.time()-t0)
    print(f"{kind:12s} a={aG:5.1f} {tname:6s} best={out['best']:6.1f} sobol_med={out['sobol_median']:6.1f} amp={out['amp']:.2e} Mmed={out['M_obj_median_kg']:.2e} n={out['n_sites']} mind={out['min_dist_from_us_gpc']:.2f} buf_viol={inside} t={out['secs']:.0f}s",flush=True)
    return out

if __name__=="__main__":
    targets=sys.argv[1].split(",") if len(sys.argv)>1 else ["alone","with4"]
    if len(sys.argv)>2 and sys.argv[2]=="bare": BUFFER=False
    tasks=[(k,aG,t) for t in targets for aG in SPACINGS for k in KINDS]
    tasks.sort(key=lambda z:z[1])     # small spacings (slow) first
    tag=f"{'_'.join(targets)}_{'buf' if BUFFER else 'bare'}"; res=[]
    with Pool(7) as p, open(f"partial_{tag}.jsonl","w") as fh:
        for o in p.imap_unordered(task,tasks,chunksize=1): res.append(o); fh.write(json.dumps(o)+"\n"); fh.flush()
    json.dump(res,open(f"search_{'_'.join(targets)}_{'buf' if BUFFER else 'bare'}.json","w"),indent=1)
