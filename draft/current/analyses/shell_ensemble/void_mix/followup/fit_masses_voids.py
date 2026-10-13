"""Fit k=1..4 (placed objects may not overlap each other) intentionally placed objects, (A) masses only (M>0) vs (B) masses or voids (sign free), same procedure.
In the ring plane, free angle and free distance 15.0-42.9 Gpc (1.05-3 R_vis, the surround region), compensated as in with_clumps
(each object's net mass removed/added back over the surround region ∝ exp-lock rho). Size cap (both signs, symmetric):
|M| <= rho_bar * (4/3) pi (s - R_vis)^3, i.e. a void may at most be completely empty (delta=-1) out to a radius that just reaches
our visible sphere; its radius grows with its deficit (R = (3|M|/(4 pi rho_bar))^(1/3) at delta=-1). Clusters get the same cap
(a delta=+1 sphere of that radius) so both signs are treated alike."""
import json, os, numpy as np
from scipy.optimize import least_squares
HERE=os.path.dirname(os.path.abspath(__file__)); WC="/workspace/clump_scan/shell_ensemble/with_clumps"
g={"__file__":os.path.join(WC,"ensemble_with_clumps.py"),"__name__":"imp"}
src=open(g["__file__"]).read().split('\nif __name__')[0]; os.chdir(WC); exec(src,g); os.chdir(HERE)
R,rms,RV,FLOOR,KG,Fc,phi=[g[k] for k in ("R","rms","RV","FLOOR","KG","Fc","phi")]
field=g["g"]["g"]["field"]

GPC=3.0857e25; RHO=2.9e-27
def kcap(D): return RHO*4/3*np.pi*(((D-1.0)*RV*GPC)**3)/KG          # cap in fitted-amplitude units
def unpack(q,k,signed):
    a=q[0:k]; D=q[k:2*k]; u=q[2*k:3*k]
    amp=kcap(D)*(np.tanh(u) if signed else 0.5*(1+np.tanh(u)))
    return a,D,amp
def model(q,k,signed):
    a,D,amp=unpack(q,k,signed); return sum(amp[i]*(field(np.radians(a[i]),D[i])-Fc) for i in range(k))
def overlaps(q,k,signed):
    a,D,amp=unpack(q,k,signed); P=np.c_[D*np.cos(np.radians(a)),D*np.sin(np.radians(a))]
    rad=(3*np.abs(amp)*KG/(4*np.pi*RHO))**(1/3)/GPC/RV          # sphere radius at |delta|=1, R_vis units
    ov=[max(0.0,rad[i]+rad[j]-np.linalg.norm(P[i]-P[j])) for i in range(k) for j in range(i+1,k)]
    return np.array(ov) if ov else np.zeros(0)
PEN=10.0   # no two placed objects may overlap (each sized at |delta|=1)
def fit(k,signed,nstart=150,seed=0):
    rng=np.random.default_rng(seed); best=None
    starts=[]
    for _ in range(nstart):
        a0=rng.uniform(0,360,k); D0=rng.uniform(1.1,2.8,k); m0=rng.uniform(1e51,3e52,k)/KG*(rng.choice([-1,1],k) if signed else 1)
        starts.append((a0,D0,m0))
    fw=json.load(open("/workspace/clump_scan/scan.json"))["sym"]["4"]["p"]            # also start from the fitted four (and subsets)
    A4=np.array([phi+fw[2],phi-fw[2],phi+fw[4],phi-fw[4]]); K4=np.array([fw[1],fw[1],fw[3],fw[3]])
    import itertools
    for sub in itertools.combinations(range(4),k): starts.append((A4[list(sub)],np.full(k,fw[0]),K4[list(sub)]))
    for a0,D0,m0 in starts:
        c=kcap(D0); x=np.clip(m0/c,-0.999,0.999) if signed else np.clip(2*m0/c-1,-0.999,0.999)
        q0=np.r_[a0,D0,np.arctanh(x)]
        try: f=least_squares(lambda q:np.r_[model(q,k,signed)-R,PEN*overlaps(q,k,signed)],q0,bounds=(np.r_[[-1e4]*k,[1.05]*k,[-8]*k],np.r_[[1e4]*k,[3.0]*k,[8]*k]),max_nfev=2000)
        except Exception: continue
        if best is None or f.cost<best.cost: best=f
    a,D,amp=unpack(best.x,k,signed); left=R-model(best.x,k,signed)
    M=amp*KG; rad=(3*np.abs(M)/(4*np.pi*RHO))**(1/3)/GPC
    return dict(rms_over_floor=rms(left)/FLOOR,ell_deg=(a%360).tolist(),dist_gpc=(D*RV).tolist(),M_kg=M.tolist(),
                kind=["void" if m<0 else "cluster" for m in M],radius_gpc_at_abs_delta_1=rad.tolist(),
                frac_of_cap=(np.abs(amp)/kcap(D)).tolist(),max_overlap_gpc=float(overlaps(best.x,k,signed).max()*RV) if k>1 else 0.0,net_kg=float(M.sum()),q=best.x.tolist())
if __name__=="__main__":
    out={}
    for k in (1,2,3,4):
        for signed in (False,True):
            r=fit(k,signed); out[f"k{k}_{'masses_voids' if signed else 'masses_only'}"]=r
            print(k,"masses+voids" if signed else "masses only",f"{r['rms_over_floor']:.1f}x",[f"{m:.2e}" for m in r["M_kg"]],[f"{d:.1f}" for d in r["dist_gpc"]],[f"{x:.0f}" for x in r["ell_deg"]],[f"{x:.2f}" for x in r["frac_of_cap"]],"overlap",f"{r['max_overlap_gpc']:.3f}",flush=True)
    json.dump(out,open("fit_masses_voids.json","w"),indent=1)
