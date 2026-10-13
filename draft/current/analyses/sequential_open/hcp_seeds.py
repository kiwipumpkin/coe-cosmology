"""Fill-seed scatter of the Delaunay link energy u and disorder entropy for PERFECT lattices with the web mix (voids in random 1/12 of
tetra holes, or 1/6 of octa holes), HCP vs FCC, 16 seeds each. Same recipe as fcc_ladder.py / hcp_lattice.py."""
import json,numpy as np
exec(open("hcp_lattice.py").read().split('if __name__=="__main__":')[0])
L0=14*A; MG=2*A; out={}
def fcc_octa(L,a,fill):
    n=int(np.ceil(L*1.8/a))+1; c,t,o=sites(0,int(L/a)+1,a); ov=o[rng.random(len(o))<fill]
    P=np.vstack([c,ov]); K=np.r_[np.zeros(len(c),int),np.ones(len(ov),int)]; m=np.all((P>=0)&(P<L),axis=1); return P[m],K[m]
for s in range(16):
    rng=np.random.default_rng(100+s)
    for nm,fn in (("HCP tetra 1/12",lambda: crystal_hcp(L0,A,1/12)),("FCC tetra 1/12",lambda: crystal_fcc(L0,A,1/12)),
                  ("HCP octa 1/6",lambda: crystal_hcp(L0,A,1/6,hole="octa")),("FCC octa 1/6",lambda: fcc_octa(L0,A,1/6))):
        P,K=fn(); r,_=run(P,K,L0,nm,MG); out.setdefault(nm,[]).append((r["u_dimless"],r["S_hist_logvol"],r["mean_degree"],r["U_cc"],r["U_cv"],r["U_vv"],r["frac_links_cv"],r["frac_links_vv"]))
S={k:dict(u_mean=float(np.mean([x[0] for x in v])),u_sd=float(np.std([x[0] for x in v])),S_hist_mean=float(np.mean([x[1] for x in v])),
          deg_mean=float(np.mean([x[2] for x in v])),frac_cv=float(np.mean([x[6] for x in v])),frac_vv=float(np.mean([x[7] for x in v])),runs=v) for k,v in out.items()}
for k,v in S.items(): print(k,{a:b for a,b in v.items() if a!="runs"})
json.dump(S,open("hcp_lattice/hcp_fill_seeds.json","w"),indent=1)
