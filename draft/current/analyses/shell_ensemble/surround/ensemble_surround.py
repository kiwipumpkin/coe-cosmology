"""Surround variant of ../ensemble.py: objects fill all space 14.3-42.9 Gpc around the OBSERVER (every direction) that is
in the CoE shell D-d..D+d or beyond D+d (out to D+3d); number density follows the exp-lock density gradient. Physics/scoring reused verbatim
from ../ensemble.py (field3, amps, summarize, LEVELS, KG, residual R from scan.py)."""
import json, os, sys, time, warnings, numpy as np
warnings.filterwarnings("ignore",category=RuntimeWarning)
SURR=os.path.dirname(os.path.abspath(__file__)); ORIG=os.path.join(os.path.dirname(SURR),"ensemble.py")
g={"__file__":ORIG,"__name__":"surround_import"}
exec(open(ORIG).read().split('if __name__=="__main__":')[0],g)   # unchanged tidal code, residual, scoring
os.chdir(SURR)
field3,amps,summarize,LEVELS,KG,R,rms,radius_gpc,RV,FLOOR,phi=[g[k] for k in
    ("field3","amps","summarize","LEVELS","KG","R","rms","radius_gpc","RV","FLOOR","phi")]

# geometry (R_vis units), same CoE/shell as original
D=58.2/RV; d=1.0                    # d = R_vis = 14.3 Gpc (CatWISE lock distance)
A_IN,B_SH=D-d,D+d                   # shell D-d..D+d = 43.9..72.5 Gpc from CoE
S_MAX=3*d                           # region: everything within 3 R_vis = 42.9 Gpc of us (2 R_vis of space around the visible sphere)
C_OUT=D+S_MAX                       # => farthest point from CoE = D+3d = 101.1 Gpc (beyond-D+d layer reaches here)
# density law: exp lock (main locked law), rho = rho0 exp(-r/L), L from delta=[rho(D-d)-rho(D+d)]/rho(D)=0.031
L_GPC=922.6487386323013; L=L_GPC/RV  # = d/asinh(0.0155); D-independent, so identical for D=48.3 or 58.2

def coe_vec(ell): return D*np.array([np.cos(np.radians(ell)),np.sin(np.radians(ell)),0.0])

def _draw(n,rng,C,s0,kind,grad):
    """n positions, number density per unit volume ∝ rho(r_CoE) (clusters, kind>0) or ∝ 1/rho (voids, kind<0)
    [∝ 1 if grad=False], over the region s0 < s <= S_MAX around the observer with r_CoE >= D-d
    (shell D-d..D+d plus everything beyond D+d; the CoE interior r<D-d stays empty as in the original)."""
    s0=np.broadcast_to(np.asarray(s0,float),(n,)).copy(); out=np.empty((n,3)); need=np.arange(n)
    while need.size:
        m=need.size*3; idx=np.repeat(need,3)
        u=rng.normal(size=(m,3)); u/=np.linalg.norm(u,axis=1)[:,None]
        s=(rng.uniform(s0[idx]**3,S_MAX**3))**(1/3); P=s[:,None]*u          # uniform in volume of the ball layer
        r=np.linalg.norm(P-C,axis=1); ok=r>=A_IN
        if grad: ok&=rng.uniform(size=m)<(np.exp(-(r-A_IN)/L) if kind>0 else np.exp((r-C_OUT)/L))
        filled=set()
        for j in np.where(ok)[0]:
            if idx[j] not in filled: out[idx[j]]=P[j]; filled.add(idx[j])
        need=np.array([i for i in need if i not in filled],int)
    return out

def realization(N,lev,rng,ell_coe,buffer,grad):
    lo,hi=LEVELS[lev]; nh=N//2
    Mc=10**rng.uniform(lo,hi,nh); Mv=10**rng.uniform(lo,hi,N-nh); Mv*=Mc.sum()/Mv.sum()   # identical to original
    M=np.r_[Mc,-Mv]
    rmin=(1.0+np.r_[radius_gpc(Mc,1.0),radius_gpc(Mv,0.5)]/RV) if buffer else np.ones(N)
    C=coe_vec(ell_coe); P=np.empty((N,3))
    P[:nh]=_draw(nh,rng,C,rmin[:nh],1,grad); P[nh:]=_draw(N-nh,rng,C,rmin[nh:],-1,grad)
    beyond=np.linalg.norm(P-C,axis=1)>B_SH
    Dd=np.linalg.norm(P,axis=1); al=np.arctan2(P[:,1],P[:,0]); be=np.arcsin(P[:,2]/Dd)
    w=(M/KG)@field3(al,be,Dd)
    return w,M,P,beyond

def run(N,lev,nreal,ell_coe,buffer,grad,seed):
    rng=np.random.default_rng(seed); A0=amps(R)
    res=dict(rms=[],rms_best=[],s_best=[],var_w=[],amps=[],amps_best=[],minD=[],frac_back=[],frac_near=[],frac_beyond=[])
    keep=None
    for i in range(nreal):
        w,M,P,bey=realization(N,lev,rng,ell_coe,buffer,grad)
        left=R-w; s=float(np.dot(R,w)/np.dot(w,w)); lb=R-s*w
        res["rms"].append(rms(left)); res["rms_best"].append(rms(lb)); res["s_best"].append(s)
        res["var_w"].append(float(np.mean(w**2))); res["amps"].append(amps(left)); res["amps_best"].append(amps(lb))
        sd=np.linalg.norm(P,axis=1); res["minD"].append(float(sd.min()*RV))
        C=coe_vec(ell_coe); res["frac_back"].append(float(np.mean(P@C<0)))     # far-side hemisphere (away from CoE)
        res["frac_near"].append(float(np.mean(sd<2.0))); res["frac_beyond"].append(float(bey.mean()))
        if i==0: keep=dict(M=M.tolist(),P=(P*RV).tolist(),beyond=bey.tolist())
    for k in res: res[k]=np.array(res[k])
    return res,keep,A0

if __name__=="__main__":
    NREAL=int(sys.argv[1]) if len(sys.argv)>1 else 300
    ELL_COE=phi+180.0
    out=dict(geometry=dict(D_gpc=58.2,d_gpc=RV,shell_gpc=[A_IN*RV,B_SH*RV],beyond_gpc=[B_SH*RV,C_OUT*RV],region_gpc_from_observer=[RV,S_MAX*RV],
             region="all space 14.3-42.9 Gpc from the observer (every direction) with CoE distance >= D-d: shell D-d..D+d plus beyond D+d",
             density_law="exp lock rho=rho0 exp(-r/L), L=%.2f Gpc (delta=0.031 over D±d)"%L_GPC,
             gradient_rule="clusters: position density ∝ rho(r_CoE); voids: ∝ 1/rho(r_CoE)"),
             base_rms_over_floor=rms(R)/FLOOR,base_amps_over_floor=(amps(R)/FLOOR).tolist(),ell_coe=ELL_COE,settings={},examples={})
    seed=1; t0=time.time()
    for grad in (True,False):
        for buffer in (False,True):
            for N in (200,500):
                for lev in "abc":
                    res,keep,A0=run(N,lev,NREAL,ELL_COE,buffer,grad,seed); seed+=1
                    key=f"{'grad' if grad else 'unif'}_{'buf' if buffer else 'bare'}_N{N}_{lev}"
                    s=summarize(res,A0,lev); s["frac_far_hemisphere"]=float(res["frac_back"].mean()); s["frac_within_28.6"]=float(res["frac_near"].mean()); s["frac_beyond"]=float(res["frac_beyond"].mean())
                    out["settings"][key]=s; out["examples"][key]=keep
                    np.savez_compressed(f"raw_{key}.npz",**res)
                    print(key,f"med={s['median']:.1f} p5={s['p5']:.1f} min={s['min']:.1f} best_med={s['best_median']:.1f} best_min={s['best_min']:.1f} "
                          f"f<237={s['frac_below_base']:.3f} f35={s['frac_35_strict']:.3f}/{s['frac_35_quad']:.3f} back={s['frac_far_hemisphere']:.2f} "
                          f"near={s['frac_within_28.6']:.3f} beyond={s['frac_beyond']:.3f} E[rms2]={s['mean_rms2_over_floor2']:.4g} pred={s['predicted_rms2_over_floor2']:.4g} minD={s['min_obj_dist_gpc']:.2f} t={time.time()-t0:.0f}s",flush=True)
    json.dump(out,open("ensemble.json","w"),indent=1)
