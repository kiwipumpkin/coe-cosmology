"""Surround variant of ../ensemble.py: objects isotropic around the OBSERVER (full 4pi), 70% in the CoE shell
D-d..D+d, 30% beyond D+d (D+d..D+3d), optional exp-lock density gradient. Physics/scoring reused verbatim
from ../ensemble.py (field3, amps, summarize, LEVELS, KG, residual R from scan.py)."""
import json, os, sys, time, warnings, numpy as np
warnings.filterwarnings("ignore",category=RuntimeWarning)
SURR=os.path.dirname(os.path.abspath(__file__)); ORIG=os.path.join(os.path.dirname(SURR),"ensemble.py")
g={"__file__":ORIG,"__name__":"surround_import"}
exec(open(ORIG).read().split('if __name__=="__main__":')[0],g)   # unchanged tidal code, residual, scoring
os.chdir(SURR)
field3,amps,summarize,LEVELS,KG,R,rms,radius_gpc,RV,FLOOR,phi=[g[k] for k in
    ("field3","amps","summarize","LEVELS","KG","R","rms","radius_gpc","RV","FLOOR","phi")]

# geometry (R_vis units), same as original
D=58.2/RV; d=1.0                    # d = R_vis = 14.3 Gpc (CatWISE lock distance)
A_IN,B_SH,C_OUT=D-d,D+d,D+3*d       # 43.9, 72.5, 101.1 Gpc: shell, and a beyond-layer of the same thickness 2d
F_BEYOND=0.30
# density law: exp lock (main locked law), rho = rho0 exp(-r/L), L from delta=[rho(D-d)-rho(D+d)]/rho(D)=0.031
L_GPC=922.6487386323013; L=L_GPC/RV  # = d/asinh(0.0155); D-independent, so identical for D=48.3 or 58.2

def coe_vec(ell): return D*np.array([np.cos(np.radians(ell)),np.sin(np.radians(ell)),0.0])

def _draw(n,rng,C,band,s0,kind,grad):
    """n positions: direction isotropic about observer; distance s along the ray restricted to CoE-radius band
    and s>s0, density along ray ∝ s^2 (uniform volume) × gradient weight (clusters ∝ rho, voids ∝ 1/rho)."""
    s0=np.broadcast_to(np.asarray(s0,float),(n,)).copy(); out=np.empty((n,3)); need=np.arange(n)
    lo_b,hi_b=band
    while need.size:
        m=need.size; u=rng.normal(size=(m,3)); u/=np.linalg.norm(u,axis=1)[:,None]
        p=u@C; S0=s0[need]
        Shi=p+np.sqrt(p*p-D*D+hi_b**2)                       # r<=hi  <=> s<=Shi
        if lo_b<D:                                           # inner edge of shell: r>=lo excludes (Sa1,Sa2)
            disc=p*p-D*D+lo_b**2; sq=np.sqrt(np.clip(disc,0,None))
            Sa1=np.where(disc>0,p-sq,np.inf); Sa2=np.where(disc>0,p+sq,-np.inf)
            iv=[(S0,np.minimum(Sa1,Shi)),(np.maximum(Sa2,S0),Shi)]
        else:                                                # beyond layer: r>=lo <=> s>=Slo
            Slo=p+np.sqrt(p*p-D*D+lo_b**2); iv=[(np.maximum(Slo,S0),Shi)]
        w=np.stack([np.clip(h**3-l**3,0,None) for l,h in iv],1); tot=w.sum(1)
        c=(rng.uniform(size=m)*tot)[:,None]>np.cumsum(w,1); k=c.sum(1)
        lo=np.choose(k,[l for l,h in iv]); hi=np.choose(k,[h for l,h in iv])
        s=(lo**3+rng.uniform(size=m)*(hi**3-lo**3))**(1/3); P=s[:,None]*u
        r=np.linalg.norm(P-C,axis=1)
        if grad:   # rejection on exp-lock density, normalised to its max inside the band
            acc=np.exp(-(r-lo_b)/L) if kind>0 else np.exp((r-hi_b)/L)
            ok=rng.uniform(size=m)<acc
        else: ok=np.ones(m,bool)
        ok&=tot>0
        out[need[ok]]=P[ok]; need=need[~ok]
    return out

def realization(N,lev,rng,ell_coe,buffer,grad):
    lo,hi=LEVELS[lev]; nh=N//2
    Mc=10**rng.uniform(lo,hi,nh); Mv=10**rng.uniform(lo,hi,N-nh); Mv*=Mc.sum()/Mv.sum()   # identical to original
    M=np.r_[Mc,-Mv]; sign=np.sign(M)
    rmin=(1.0+np.r_[radius_gpc(Mc,1.0),radius_gpc(Mv,0.5)]/RV) if buffer else np.ones(N)
    # 30% of clusters and 30% of voids go beyond D+d (random members)
    beyond=np.zeros(N,bool)
    for idx in (np.arange(nh),np.arange(nh,N)):
        beyond[rng.choice(idx,int(round(F_BEYOND*idx.size)),replace=False)]=True
    C=coe_vec(ell_coe); P=np.empty((N,3))
    for bnd,sel in (((A_IN,B_SH),~beyond),((B_SH,C_OUT),beyond)):
        for kind in (1,-1):
            ii=np.where(sel&(sign==kind))[0]
            if ii.size: P[ii]=_draw(ii.size,rng,C,bnd,rmin[ii],kind,grad)
    Dd=np.linalg.norm(P,axis=1); al=np.arctan2(P[:,1],P[:,0]); be=np.arcsin(P[:,2]/Dd)
    w=(M/KG)@field3(al,be,Dd)
    return w,M,P,beyond

def run(N,lev,nreal,ell_coe,buffer,grad,seed):
    rng=np.random.default_rng(seed); A0=amps(R)
    res=dict(rms=[],rms_best=[],s_best=[],var_w=[],amps=[],amps_best=[],minD=[],frac_back=[],frac_near=[])
    keep=None
    for i in range(nreal):
        w,M,P,bey=realization(N,lev,rng,ell_coe,buffer,grad)
        left=R-w; s=float(np.dot(R,w)/np.dot(w,w)); lb=R-s*w
        res["rms"].append(rms(left)); res["rms_best"].append(rms(lb)); res["s_best"].append(s)
        res["var_w"].append(float(np.mean(w**2))); res["amps"].append(amps(left)); res["amps_best"].append(amps(lb))
        sd=np.linalg.norm(P,axis=1); res["minD"].append(float(sd.min()*RV))
        C=coe_vec(ell_coe); res["frac_back"].append(float(np.mean(P@C<0)))     # far-side hemisphere (away from CoE)
        res["frac_near"].append(float(np.mean(sd<2.0)))                        # within 28.6 Gpc of observer
        if i==0: keep=dict(M=M.tolist(),P=(P*RV).tolist(),beyond=bey.tolist())
    for k in res: res[k]=np.array(res[k])
    return res,keep,A0

if __name__=="__main__":
    NREAL=int(sys.argv[1]) if len(sys.argv)>1 else 300
    ELL_COE=phi+180.0
    out=dict(geometry=dict(D_gpc=58.2,d_gpc=RV,shell_gpc=[A_IN*RV,B_SH*RV],beyond_gpc=[B_SH*RV,C_OUT*RV],f_beyond=F_BEYOND,
             directions="isotropic about observer (4pi)",radial="uniform in volume along each ray (p ∝ s^2) within the band, s > R_vis",
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
                    s=summarize(res,A0,lev); s["frac_far_hemisphere"]=float(res["frac_back"].mean()); s["frac_within_28.6"]=float(res["frac_near"].mean())
                    out["settings"][key]=s; out["examples"][key]=keep
                    np.savez_compressed(f"raw_{key}.npz",**res)
                    print(key,f"med={s['median']:.1f} p5={s['p5']:.1f} min={s['min']:.1f} best_med={s['best_median']:.1f} best_min={s['best_min']:.1f} "
                          f"f<237={s['frac_below_base']:.3f} f35={s['frac_35_strict']:.3f}/{s['frac_35_quad']:.3f} back={s['frac_far_hemisphere']:.2f} "
                          f"near={s['frac_within_28.6']:.3f} E[rms2]={s['mean_rms2_over_floor2']:.4g} pred={s['predicted_rms2_over_floor2']:.4g} minD={s['min_obj_dist_gpc']:.2f} t={time.time()-t0:.0f}s",flush=True)
    json.dump(out,open("ensemble.json","w"),indent=1)
