"""Translate each 1x scenario's required excess mass (M+) and deficit (M-) into numbers of realistic objects; feasibility numbers."""
import json, numpy as np
MSUN=1.989e30; GPC=3.0857e25; MPC=GPC/1e3; RHO=2.9e-27
rho_mpc=RHO*MPC**3                                   # kg per Mpc^3 at mean density (8.5e40 kg = 4.3e10 Msun)
def sph(R): return 4/3*np.pi*R**3
CL=dict(name="rich cluster",excess=1e15*MSUN,R=2.0,delta=200)
SCL=dict(name="supercluster (Shapley-like)",R=50.0,delta=1.2); SCL["excess"]=SCL["delta"]*rho_mpc*sph(50)
VO=dict(name="void",R=28.6,delta=-0.4); VO["deficit"]=0.4*rho_mpc*sph(28.6)
SVO=dict(name="supervoid (Eridanus-like)",R=200.0,delta=-0.2); SVO["deficit"]=0.2*rho_mpc*sph(200)
# natural content (local_web/volumes_counts.md): excess mass in clusters / deficit in voids, and counts
NAT={"inside":dict(V=1.22e4,clu_mass=3.0e53,void_def=5.7e53,n_clu=6.6e8,n_void=9.8e7),
     "outside":dict(V=8.82e4,clu_mass=2.2e54,void_def=4.1e54,n_clu=4.8e9,n_void=7.1e8)}   # outside = tidal region 14.3-30 Gpc
NAT["both"]={k:NAT["inside"][k]+NAT["outside"][k] for k in NAT["inside"]}
F=json.load(open("fields.json"))
try: IA=json.load(open("inside_asym_best.json"))
except Exception: IA=None
S=[]
def add(geom,sym,kind,label,Mp,Mm,left,extra):
    S.append(dict(geom=geom,sym=sym,kind=kind,label=label,M_plus=Mp,M_minus=Mm,left=left,
        n_clusters=Mp/CL["excess"],n_superclusters=Mp/SCL["excess"],n_voids=Mm/VO["deficit"],n_supervoids=Mm/SVO["deficit"],
        nat_cluster_equiv=NAT[geom]["clu_mass"]/CL["excess"],nat_void_equiv=NAT[geom]["void_def"]/VO["deficit"],
        frac_of_natural_cluster_mass=Mp/NAT[geom]["clu_mass"],frac_of_natural_void_deficit=Mm/NAT[geom]["void_def"],**extra))
# smooth fields (2 Gpc smoothing, min-norm)
for key,geom,sym,lab in [("inside_mirror","inside","mirror","inside, smooth field"),("outside_mirror","outside","mirror","outside, smooth field"),
                         ("both_mirror","both","mirror","both, smooth field"),("both_staggered_in+out-","both","no-mirror","both, staggered field (dense side inside, empty side outside)"),
                         ("inside_onesided+","inside","no-mirror","inside, one-sided field"),("inside_staggered_near+far-","inside","no-mirror","inside, staggered field (dense side near ring, empty side nearer us)"),("outside_onesided+","outside","no-mirror","outside, one-sided field")]:
    f=F[key]
    if not f["reach"]: S.append(dict(geom=geom,sym=sym,kind="field",label=lab,reach=False)); continue
    p=f["parts"]; ex=dict(reach=True,peak_delta=max(q["peak"] for q in p),rms_delta=[q["rms"] for q in p],mean_delta_over=[q["mean_d_over"] for q in p],mean_delta_under=[q["mean_d_under"] for q in p])
    ex['dmin']=min(q['dmin'] for q in p); ex['dmax']=max(q['dmax'] for q in p)
    add(geom,sym,"field",lab,f["M_plus"],f["M_minus"],f["left"],ex)
# point-mass fits (net-zero mixed bests)
add("inside","mirror","points","inside, 4 masses + 4 voids (0.18x)",5.4e52,5.4e52,0.18,dict(reach=True,core_delta=100,void_delta=-0.7,void_R_gpc="3.5-4"))
add("outside","mirror","points","outside, 6 masses + 5 voids (0.58x)",6.77e52,6.76e52,0.58,dict(reach=True,core_delta=100,void_delta=-0.7,void_R_gpc="3.5-4"))
add("outside","no-mirror","points","outside, 7 masses + 5 voids, no mirror (0.73x)",5.82e52,5.81e52,0.73,dict(reach=True,core_delta=100,void_delta=-0.7,void_R_gpc="2.3-4"))
if IA: add("inside","no-mirror","points",IA["label"],IA["M_plus"],IA["M_minus"],IA["left"],dict(reach=True))
# feasibility: point object of 2e52 kg spread into a lobe of radius R (Gpc) -> mean delta; deep void Rv=4 Gpc delta=-0.7 filled with real voids
feas={}
for Rl in (0.8,1,2,3):
    feas[f"2e52kg_in_R{Rl}Gpc_delta"]=2e52/(RHO*sph(Rl)*GPC**3)
feas["void_fill_needed_for_delta-0.7_with_delta-0.4_voids"]=0.7/0.4
feas["void_fill_needed_for_delta-0.7_with_supervoids_-0.2"]=0.7/0.2
consts=dict(rich_cluster_kg=CL["excess"],supercluster_excess_kg=SCL["excess"],void_deficit_kg=VO["deficit"],supervoid_deficit_kg=SVO["deficit"],rho_kg_per_mpc3=rho_mpc)
json.dump(dict(scenarios=S,feasibility=feas,constants=consts,natural=NAT),open("counts.json","w"),indent=1)
for s in S:
    if not s.get("reach",True) or "M_plus" not in s: print(f"{s['label']:70s} cannot reach 1x"); continue
    print(f"{s['label']:70s} M+={s['M_plus']:.1e} clusters={s['n_clusters']:.1e} SC={s['n_superclusters']:.1e} voids={s['n_voids']:.1e} SV={s['n_supervoids']:.1e} "
          f"of natural: clu {s['frac_of_natural_cluster_mass']*100:.2g}% void {s['frac_of_natural_void_deficit']*100:.2g}%  peak={s.get('peak_delta','-')}")
print(feas); print(consts)
