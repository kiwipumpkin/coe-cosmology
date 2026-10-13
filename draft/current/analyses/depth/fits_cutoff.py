"""(A) Best leftover vs R_cut: refit k placed objects with distances restricted to [15.0 Gpc, R_cut].
Same rules as ../shell_ensemble/void_mix/followup/fit_masses_voids.py (imported verbatim: model, compensation, size cap,
no-overlap penalty); only the distance upper bound changes. Optional: keep objects >= D-d (43.9 Gpc) from the CoE.
usage: python3 fits_cutoff.py <k> <signed 0/1> <coe 0/1> <tag>"""
import json, os, sys, numpy as np, itertools
from scipy.optimize import least_squares
HERE=os.path.dirname(os.path.abspath(__file__)); FU="/workspace/clump_scan/shell_ensemble/void_mix/followup"
fm={"__file__":os.path.join(FU,"fit_masses_voids.py"),"__name__":"imp"}
exec(open(fm["__file__"]).read().split('\nif __name__')[0],fm); os.chdir(HERE)
R,rms,RV,FLOOR,KG,phi,model,overlaps,kcap,unpack,PEN=[fm[x] for x in ("R","rms","RV","FLOOR","KG","phi","model","overlaps","kcap","unpack","PEN")]
C2=fm["g"]["C"][:2]; A_IN=fm["g"]["g"]["A_IN"]
RCUTS=[16,18,20,22.5,25,28.6,32,36,42.9,50,57.2]
DMIN=1.05                                            # 15.0 Gpc, as in the earlier fits (cap -> 0 at 14.3 anyway)
def coe_pen(q,k,signed):
    a,D,amp=unpack(q,k,signed); P=np.c_[D*np.cos(np.radians(a)),D*np.sin(np.radians(a))]
    return np.maximum(0.0,A_IN-np.linalg.norm(P-C2,axis=1))
def fit(k,signed,coe,Dmax,warm,nstart=120,seed=0):
    rng=np.random.default_rng(seed); starts=[]
    for _ in range(nstart):
        a0=rng.uniform(0,360,k); D0=rng.uniform(DMIN+0.01,max(DMIN+0.02,min(2.9,Dmax)-0.01),k)
        m0=rng.uniform(1e51,3e52,k)/KG*(rng.choice([-1,1],k) if signed else 1)
        c=kcap(D0); x=np.clip(m0/c,-0.999,0.999) if signed else np.clip(2*m0/c-1,-0.999,0.999)
        starts.append(np.r_[a0,D0,np.arctanh(x)])
    for w in warm:                                   # warm starts: previous R_cut best + unrestricted fit, distances clipped
        w=np.array(w,float).copy(); w[k:2*k]=np.clip(w[k:2*k],DMIN,Dmax); starts.append(w)
    lb=np.r_[[-1e4]*k,[DMIN]*k,[-8]*k]; ub=np.r_[[1e4]*k,[Dmax]*k,[8]*k]
    def fun(q):
        r=[model(q,k,signed)-R,PEN*overlaps(q,k,signed)]
        if coe: r.append(PEN*coe_pen(q,k,signed))
        return np.concatenate(r)
    best=None
    for q0 in starts:
        q0=np.clip(q0,lb+1e-9,ub-1e-9)
        try: f=least_squares(fun,q0,bounds=(lb,ub),max_nfev=2000)
        except Exception: continue
        if best is None or f.cost<best.cost: best=f
    q=best.x; a,D,amp=unpack(q,k,signed); left=R-model(q,k,signed)
    P=np.c_[D*np.cos(np.radians(a)),D*np.sin(np.radians(a))]
    return dict(R_cut_gpc=Dmax*RV,rms_over_floor=rms(left)/FLOOR,M_kg=(amp*KG).tolist(),dist_gpc=(D*RV).tolist(),ell_deg=(a%360).tolist(),
                dist_from_CoE_gpc=(np.linalg.norm(P-C2,axis=1)*RV).tolist(),frac_of_cap=(np.abs(amp)/kcap(D)).tolist(),
                max_overlap_gpc=float(overlaps(q,k,signed).max()*RV) if k>1 else 0.0,
                max_coe_violation_gpc=float(coe_pen(q,k,signed).max()*RV),q=q.tolist())
if __name__=="__main__":
    k=int(sys.argv[1]); signed=bool(int(sys.argv[2])); coe=bool(int(sys.argv[3])); tag=sys.argv[4]
    prev=json.load(open(os.path.join(FU,"fit_masses_voids.json")))[f"k{k}_{'masses_voids' if signed else 'masses_only'}"]["q"]
    out=[]; warm=[prev]
    for Rc in RCUTS:
        r=fit(k,signed,coe,Rc/RV,warm); out.append(r); warm=[prev,r["q"]]
        print(tag,f"R_cut={Rc:5.1f}  leftover={r['rms_over_floor']:.1f}x  dist={np.round(r['dist_gpc'],1).tolist()}  M={[f'{m:.1e}' for m in r['M_kg']]}  rCoE={np.round(r['dist_from_CoE_gpc'],1).tolist()}",flush=True)
        json.dump(dict(k=k,signed=signed,coe_constraint=coe,scan=out),open(f"fits_{tag}.json","w"),indent=1)
