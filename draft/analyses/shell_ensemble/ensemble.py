"""Random cluster+void web in a shell around the CoE; same tidal code as the four-clump fit (scan.py)."""
import json, os, sys, time, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); CS=os.path.dirname(HERE); os.chdir(CS)
exec(open("scan.py").read().split("out={")[0])          # R, th, phi, ELL, field, modes, rms  (unchanged)
os.chdir(HERE)
KG=7.7e51/0.0039          # kg per unit fitted amplitude k, exactly as make_figs.py / pos_free.py / web_floor.py
FLOOR=1e-5; RV=14.3       # Gpc
GPC=3.0857e25; RHO=2.9e-27  # m, kg/m^3 (Omega_m 0.31, H0 70)
def radius_gpc(M,delta): return (3*abs(M)/(4*np.pi*RHO*abs(delta)))**(1/3)/GPC

def field3(alpha,beta,D):
    """3D generalisation of field(): unit mass at distance D (R_vis units), longitude alpha, latitude beta
    (rad) out of the ring plane; tidal radial acceleration on the unit ring of in-plane sightlines,
    minus the uniform (observer) term, minus the mean. Arrays alpha,beta,D -> (n_obj, n_theta)."""
    alpha,beta,D=np.atleast_1d(alpha)[:,None],np.atleast_1d(beta)[:,None],np.atleast_1d(D)[:,None]
    x,y=np.cos(th)[None],np.sin(th)[None]
    X,Y,Z=D*np.cos(beta)*np.cos(alpha),D*np.cos(beta)*np.sin(alpha),D*np.sin(beta)
    dx,dy,dz=X-x,Y-y,Z; r3=(dx*dx+dy*dy+dz*dz)**1.5
    ar=(dx*x+dy*y)/r3-(X*x+Y*y)/D**3
    return ar-ar.mean(axis=1,keepdims=True)

def amps(v):  # m=1..12 amplitudes, same definition as scan.py modes()/fourier()
    m=np.arange(1,13)[:,None]
    return np.hypot(2*np.mean(v[...,None,:]*np.cos(m*th),-1),2*np.mean(v[...,None,:]*np.sin(m*th),-1))

# ---------- beta=0 check
rng=np.random.default_rng(0); chk=0.0
for _ in range(200):
    a,D=rng.uniform(0,2*np.pi),rng.uniform(1.05,10)
    chk=max(chk,np.max(np.abs(field3(a,0.0,D)[0]-field(a,D)))/np.max(np.abs(field(a,D))))
print("beta=0 check: max relative |field3-field| =",chk)

# ---------- geometry (code convention: theta=ell-phi=0 points AWAY from O, so O lies at ell=phi+180)
D_COE=58.2/RV; R_IN,R_OUT=43.9/RV,72.5/RV
def coe_vec(ell_coe): return D_COE*np.array([np.cos(np.radians(ell_coe)),np.sin(np.radians(ell_coe)),0.0])

def sample(n,rng,ell_coe,rmin_obs):
    """n points uniform in shell volume around CoE, with distance from observer > rmin_obs (array or scalar)."""
    C=coe_vec(ell_coe); out=np.empty((0,3)); rmin=np.broadcast_to(np.asarray(rmin_obs,float),(n,)).copy(); need=np.arange(n); P=np.empty((n,3))
    while need.size:
        k=need.size*2
        r=(rng.uniform(R_IN**3,R_OUT**3,k))**(1/3); u=rng.normal(size=(k,3)); u/=np.linalg.norm(u,axis=1)[:,None]
        p=C+r[:,None]*u; ok=np.linalg.norm(p,axis=1)
        take=[]; j=0
        for i in range(k):
            if j>=need.size: break
            if ok[i]>rmin[need[j]]: P[need[j]]=p[i]; j+=1
        need=need[j:]
    return P

LEVELS={"a":(48.0,49.0),"b":(49.5,50.5),"c":(51.0,52.0)}   # log10 |M| (kg), log-uniform per object
def realization(N,lev,rng,ell_coe,buffer):
    lo,hi=LEVELS[lev]; nh=N//2
    Mc=10**rng.uniform(lo,hi,nh); Mv=10**rng.uniform(lo,hi,N-nh); Mv*=Mc.sum()/Mv.sum()  # zero net contrast
    M=np.r_[Mc,-Mv]
    if buffer:  # whole object outside the vision sphere: centre > R_vis + own radius (delta=+1 clusters, -0.5 voids)
        rad=np.r_[radius_gpc(Mc,1.0),radius_gpc(Mv,0.5)]/RV; rmin=1.0+rad
    else: rmin=1.0
    P=sample(N,rng,ell_coe,rmin)
    D=np.linalg.norm(P,axis=1); al=np.arctan2(P[:,1],P[:,0]); be=np.arcsin(P[:,2]/D)
    w=(M/KG)@field3(al,be,D)          # model pattern in Delta T/T (k = M/KG per object)
    return w,M,P

def run(N,lev,nreal,ell_coe,buffer,seed):
    rng=np.random.default_rng(seed); A0=amps(R)
    res=dict(rms=[],rms_best=[],s_best=[],var_w=[],amps=[],amps_best=[],minD=[])
    keep=None
    for i in range(nreal):
        w,M,P=realization(N,lev,rng,ell_coe,buffer)
        left=R-w                                   # scan.py convention: residual = R - model
        s=float(np.dot(R,w)/np.dot(w,w)); lb=R-s*w   # best global multiplier (least squares)
        res["rms"].append(rms(left)); res["rms_best"].append(rms(lb)); res["s_best"].append(s)
        res["var_w"].append(float(np.mean(w**2))); res["amps"].append(amps(left)); res["amps_best"].append(amps(lb))
        res["minD"].append(float(np.min(np.linalg.norm(P,axis=1))*RV))
        if i==0: keep=dict(M=M.tolist(),P=(P*RV).tolist())
    for k in res: res[k]=np.array(res[k])
    return res,keep,A0

def summarize(res,A0,lev):
    r=res["rms"]/FLOOR; rb=res["rms_best"]/FLOOR; base=rms(R)/FLOOR
    def frac35(A):
        strict=np.all(A[:,2:5]<A0[2:5],1)&np.all(A[:,0:2]<=A0[0:2],1)
        quad=(np.linalg.norm(A[:,2:5],axis=1)<np.linalg.norm(A0[2:5]))&(np.linalg.norm(A[:,0:2],axis=1)<=np.linalg.norm(A0[0:2]))
        return float(strict.mean()),float(quad.mean())
    lo,hi=LEVELS[lev]
    return dict(n=len(r),median=float(np.median(r)),p5=float(np.percentile(r,5)),min=float(r.min()),
        best_median=float(np.median(rb)),best_p5=float(np.percentile(rb,5)),best_min=float(rb.min()),
        frac_below_base=float(np.mean(r<base)),frac_best_below_base=float(np.mean(rb<base)),
        frac_35_strict=frac35(res["amps"])[0],frac_35_quad=frac35(res["amps"])[1],
        frac_35_best_strict=frac35(res["amps_best"])[0],frac_35_best_quad=frac35(res["amps_best"])[1],
        mean_rms2_over_floor2=float(np.mean(res["rms"]**2))/FLOOR**2,
        predicted_rms2_over_floor2=float(rms(R)**2+np.mean(res["var_w"]))/FLOOR**2,
        median_web_rms_over_floor=float(np.median(np.sqrt(res["var_w"])))/FLOOR,
        median_abs_s_best=float(np.median(np.abs(res["s_best"]))),
        median_implied_mass_kg=float(np.median(np.abs(res["s_best"])))*10**((lo+hi)/2),
        median_amps_over_floor=(np.median(res["amps"],0)/FLOOR).tolist(),
        median_amps_best_over_floor=(np.median(res["amps_best"],0)/FLOOR).tolist(),
        min_obj_dist_gpc=float(res["minD"].min()))

if __name__=="__main__":
    NREAL=int(sys.argv[1]) if len(sys.argv)>1 else 300
    ELL_COE=phi+180.0
    out=dict(beta0_check=chk,KG_per_k=KG,base_rms_over_floor=rms(R)/FLOOR,base_amps_over_floor=(amps(R)/FLOOR).tolist(),
             ell_coe=ELL_COE,phi=phi,settings={},examples={})
    out["radii"]={lev:{"M_lo":10**lo,"M_hi":10**hi,"R_cluster_d+1_gpc":[radius_gpc(10**lo,1),radius_gpc(10**hi,1)],
                       "R_void_d-0.5_gpc":[radius_gpc(10**lo,.5),radius_gpc(10**hi,.5)]} for lev,(lo,hi) in LEVELS.items()}
    out["radius_check_7.6e51"]=[radius_gpc(7.6e51,1),radius_gpc(7.6e51,.5)]
    seed=1; t0=time.time()
    for buffer in (False,True):
        for N in (200,500):
            for lev in "abc":
                res,keep,A0=run(N,lev,NREAL,ELL_COE,buffer,seed); seed+=1
                key=f"{'buf' if buffer else 'bare'}_N{N}_{lev}"
                out["settings"][key]=summarize(res,A0,lev); out["examples"][key]=keep
                np.savez_compressed(f"raw_{key}.npz",**res)
                s=out["settings"][key]
                print(key,f"med={s['median']:.1f} p5={s['p5']:.1f} best_med={s['best_median']:.1f} best_p5={s['best_p5']:.1f} "
                      f"f<237={s['frac_below_base']:.3f} f35={s['frac_35_strict']:.3f}/{s['frac_35_quad']:.3f} "
                      f"E[rms2]={s['mean_rms2_over_floor2']:.4g} pred={s['predicted_rms2_over_floor2']:.4g} minD={s['min_obj_dist_gpc']:.2f} t={time.time()-t0:.0f}s",flush=True)
    # direction check: CoE placed at ell=phi (15.2 deg) instead
    res,keep,A0=run(500,"c",NREAL,phi,False,999); out["settings"]["bare_N500_c_CoE_at_15.2"]=summarize(res,A0,"c")
    s=out["settings"]["bare_N500_c_CoE_at_15.2"]; print("CoE@15.2 bare_N500_c",f"med={s['median']:.1f} p5={s['p5']:.1f} best_med={s['best_median']:.1f}")
    json.dump(out,open("ensemble.json","w"),indent=1)
