"""Compensated random starting web in the real 3D modeled region (Chris, 10 Oct).
Region (SEQ_OPEN=tidalcore, same as the minedit_ent runs): ball of 30 Gpc around us incl. the part of the CoE core inside it;
whole spheres inside 30 Gpc, >= 0.5 Gpc from us, no sphere straddles the CMB ring (wholly inside 14.2 or outside 14.4 Gpc).
Distributions:
  void radii  : SHAPE of the project's void size function (lss.py: Jennings, Li & Hu 2013 Vdn on Sheth & van de Weygaert 2004,
                same sigma(R) as density_vs_scale), sampled above 10 Mpc, scaled so the median is 1.5 Gpc (the body scale of the
                foam webs). Range ~1.35-2.5 Gpc.
  void depth  : Clampitt & Jain 2015 (arXiv:1404.1834) shallow voids: mean contrast inside Rv between -0.25 (edge) and -0.45; uniform
                in [-0.45,-0.25] (mean -0.35, same as shallow_core.py).
  deficit     : M_def = |delta_v| rho(x) (4/3) pi R^3 (rho = project exp-lock density at the void centre)
  cluster mass: SHAPE of the Tinker et al. 2008 mass function above 3e15 Msun (abundance matched, lss/vmn_gen, M/Mmin ~1.0-1.7),
                scaled so the median is 1.40e50 kg (the foam webs' median cluster). Excess-mass convention.
  cluster size: overdensity delta_c log-uniform 30-300 (median ~95; the foam webs' 0.158 Gpc cores are delta ~100):
                r_c = (3 M_c / (4 pi delta_c rho))^(1/3)
  compensation: clusters are added to each void until they hold M_def (last one trimmed) -> N ~ 0.9 + M_def / <M_c>;
                placed uniform-in-volume in the wall R..1.3R around their void (Hamaus 2014 / Clampitt-Jain ridge), random directions.
  centres     : uniform in volume, random sequential addition; void sizes drawn at random from the size function (enough for the target
                volume +15% spare), placed largest first; a void is kept only if all its clusters fit
                (no overlaps with any sphere, gap 0.02 Gpc); stop at void volume fraction 30% (as in the foam webs) or jam.
Writes final_comp_s<seed>.npz (same keys as the growth runs + grp = owning void of each body)."""
import sys,json,numpy as np
from scipy.spatial import cKDTree
from common import *
seed=int(sys.argv[1]); rng=np.random.default_rng(1000+seed)
SH=np.load(f"{OUT}/compweb_shapes.npz")
sys.path.insert(0,"/workspace/clump_scan/void_mass_network"); import lss
RT=30.0; GAPU=0.5; GR=0.1; GAP=0.02; VF=0.30
def vol(r): return 4/3*np.pi*r**3
rv_shape=lss.sample_void_radii(np.random.default_rng(seed),200000,10.0)/10.0; RSC=1.5/np.median(rv_shape)
mc_shape=SH["MC"]/3e15; MSC=1.40e50/np.median(mc_shape)
def ok_region(x,r):
    s=np.linalg.norm(x); return (s+r<=RT) and (s-r>=GAPU) and ((s+r<=RMIN-GR) or (s-r>=RMIN+GR))
class Hash:
    def __init__(s,h=3.0): s.h=h; s.g={}; s.X=[]; s.r=[]; s.rmax=0.0
    def free(s,x,r):
        n=int(np.ceil((r+s.rmax+GAP)/s.h)); k=np.floor(x/s.h).astype(int)
        for i in range(k[0]-n,k[0]+n+1):
            for j in range(k[1]-n,k[1]+n+1):
                for l in range(k[2]-n,k[2]+n+1):
                    for q in s.g.get((i,j,l),()):
                        d=s.X[q]-x
                        if d@d<(r+s.r[q]+GAP)**2: return False
        return True
    def add(s,x,r):
        s.g.setdefault(tuple(np.floor(x/s.h).astype(int)),[]).append(len(s.X)); s.X.append(np.asarray(x,float)); s.r.append(float(r)); s.rmax=max(s.rmax,r)
VREG=vol(RT); Rs=rv_shape[rng.choice(len(rv_shape),20000,replace=False)]*RSC; k=int(np.searchsorted(np.cumsum(vol(Rs)),1.15*VF*VREG))+1; Rlist=np.sort(Rs[:k])[::-1]   # random draw from the size function (15% spare), placed largest first
H=Hash(); X=[];r=[];M=[];kind=[];grp=[];rows=[]; vvol=0.0; g=0; fails=0
for Rv_ in Rlist:
    R=Rv_
    if vvol>=VF*VREG or fails>4000: break
    dv=rng.uniform(-0.45,-0.25); placed=False
    for t in range(60):
        c=rng.normal(size=3); c/=np.linalg.norm(c); c*=RT*rng.uniform()**(1/3)
        if not ok_region(c,R) or not H.free(c,R): continue
        Md=abs(dv)*float(seq_core.rho_gpc(c[None])[0])*vol(R)*seq_core.GPC**3; rho=float(seq_core.rho_gpc(c[None])[0])
        cl=[]; acc=0.0; Hl=[]; bad=False; tries=0
        while acc<Md*(1-1e-9):
            mc=float(mc_shape[rng.integers(len(mc_shape))]*MSC); dc=float(np.exp(rng.uniform(np.log(30),np.log(300))))
            if acc+mc>Md: mc=Md-acc
            rc=(3*mc/(4*np.pi*dc*rho*seq_core.GPC**3))**(1/3); ok=False
            for t2 in range(80):
                u=rng.normal(size=3); u/=np.linalg.norm(u); rr=(rng.uniform((R+rc+GAP)**3,(1.3*R)**3))**(1/3); x=c+u*rr
                if not ok_region(x,rc) or not H.free(x,rc): continue
                if any(np.linalg.norm(x-y)<rc+ry+GAP for y,ry in Hl): continue
                ok=True; break
            if not ok: bad=True; break
            cl.append((x,rc,mc,dc)); Hl.append((x,rc)); acc+=mc
        if bad: continue
        H.add(c,R); X.append(c); r.append(R); M.append(-Md); kind.append("void"); grp.append(g)
        for x,rc,mc,dc in cl: H.add(x,rc); X.append(x); r.append(rc); M.append(mc); kind.append("cluster"); grp.append(g)
        rows.append(dict(R=float(R),delta_v=float(dv),M_def=Md,n_clusters=len(cl),mean_Mc=acc/len(cl),delta_c=[c_[3] for c_ in cl]))
        vvol+=vol(R); g+=1; placed=True; break
    fails=0 if placed else fails+1
    if g%200==0 and placed: print(g,"voids",len(X),"bodies vfrac %.3f"%(vvol/VREG),flush=True)
X=np.array(X);r=np.array(r);M=np.array(M);kind=np.array(kind);grp=np.array(grp)
w=np.zeros(360)
for i in range(0,len(X),3000): w+=M[i:i+3000]@seq_core.unit_pattern_fast(X[i:i+3000])
from common import R as RES_CMB
left=rms(RES_CMB-w)/FLOOR
np.savez(f"{OUT}/final_comp_s{seed}.npz",Xp=X,rp=r,kp=kind,Mp=M,w=w,left=left,X=X,M=M,kind=kind,grp=grp)
Nc=np.array([q["n_clusters"] for q in rows]); Md=np.array([q["M_def"] for q in rows])
json.dump(dict(seed=seed,n_voids=int(g),n_clusters=int((kind=="cluster").sum()),void_vol_frac=vvol/VREG,left=float(left),net_mass=float(M.sum()),
     R_pct=np.percentile(r[kind=="void"],[5,50,95]).tolist(),Mc_pct=np.percentile(M[kind=="cluster"],[5,50,95]).tolist(),rc_pct=np.percentile(r[kind=="cluster"],[5,50,95]).tolist(),
     clusters_per_void_mean=float(Nc.mean()),fit_linear=np.polyfit(Md,Nc,1).tolist(),RSC=RSC,MSC=MSC,rows=rows),open(f"{OUT}/compweb_s{seed}_gen.json","w"),indent=1)
print("DONE seed",seed,"voids",g,"clusters",(kind=="cluster").sum(),"vfrac %.3f"%(vvol/VREG),"left %.1f"%left,"net %.2e"%M.sum())
