import sys,json,time; sys.path.insert(0,".")
from tri2d_core import *
from scipy.ndimage import gaussian_filter
G=6.674e-11; GPC=3.0857e25; MC=1.40e50; MV=-8.39e50; MVr=MV/MC
Lx,Ly=44,46; X,I,J,NB,E,W=lattice(Lx,Ly); N=len(X); rng=np.random.default_rng(3)
t0=time.time(); P=phi_matrix(X,I,J,W,Lx,Ly); print("phi",time.time()-t0,"s")
half=N//2
def fix(s,score):   # exact 50/50: the half with the highest score are clusters
    o=np.argsort(-score); c=np.zeros(N,bool); c[o[:half]]=True; return c
def clumps(sig,seed):
    g=np.random.default_rng(seed).normal(size=(Ly,Lx)); g=gaussian_filter(g,sig,mode="wrap"); return fix(None,g.ravel()+1e-9*np.random.default_rng(seed).normal(size=N))
CONF={}
CONF["random 50/50"]=fix(None,rng.random(N))
CONF["alternating rows (= best alternation)"]=(J%2==0)
CONF["double-row stripes"]=((J//2)%2==0)
CONF["small clumps"]=clumps(1.2,5)
CONF["big clumps"]=clumps(4.0,6)
CONF["phase separation, row interface"]=(J<Ly//2)
CONF["phase separation, column interface"]=(I<Lx//2)
def masses(c,case): return np.where(c,1.0,-1.0) if case=="A" else np.where(c,1.0,MVr)
# Monte Carlo (nonlocal Kawasaki swaps of one cluster and one void at fixed composition), anneal from random
Tn=np.geomspace(3.0,0.01,60); Tf=np.geomspace(5.0,0.01,40)
for case in ("A","B"):
    sc=1.0 if case=="A" else ((1-MVr)/2)**2
    m=masses(CONF["random 50/50"],case).copy(); m=kawasaki_nn(m,NB,Tn*sc,40*N,11); CONF[f"MC anneal, nearest-neighbour links [{case}]"]=(m>0)
    m=masses(CONF["random 50/50"],case).copy(); t0=time.time(); m=kawasaki_full(m,P,Tf*sc,15*N,12); print("full MC",case,time.time()-t0); CONF[f"MC anneal, all pairs (Ewald) [{case}]"]=(m>0)
def stats(c,case):
    m=masses(c,case); like=(c[E[:,0]]==c[E[:,1]]); cc=c[E[:,0]]&c[E[:,1]]; vv=~c[E[:,0]]&~c[E[:,1]]
    return dict(n_clusters=int(c.sum()),like_links_per_node=float(2*like.sum()/N),unlike_links_per_node=float(2*(~like).sum()/N),
                cc_per_node=float(2*cc.sum()/N),vv_per_node=float(2*vv.sum()/N),U_nn_per_node=float(U_nn(m,E)/N),U_all_per_node=float(U_full(m,P)/N))
RES={"_setup":dict(lattice="triangular, periodic box %dx%d sites = %d (box %.1f a x %.1f a)"%(Lx,Ly,N,W[0],W[1]),
     units="energies in G*M_c^2/a per node; M_c=1.40e50 kg; case A |M_v|=M_c; case B M_v=-8.39e50 kg (= %.3f M_c)"%MVr,
     all_pairs="2D-periodic Ewald sum of 1/r with neutralizing uniform background (k=0 term); includes self-image terms (constant at fixed composition)",
     MC="nonlocal Kawasaki swaps (random cluster <-> random void), geometric anneal; NN: T 3->0.01 (x scale), 60 steps x 40N attempts; all pairs: T 5->0.01, 40 x 15N attempts; B temperatures scaled by ((1-Mv/Mc)/2)^2=%.2f"%(((1-MVr)/2)**2),
     triangle_area_a2=float(np.sqrt(3)/4),link_length_a=1.0,connections_per_node=6)}
for case in ("A","B"):
    for k,c in CONF.items():
        if k.startswith("MC") and not k.endswith(f"[{case}]"): continue
        RES[f"{case} | {k}"]=stats(c,case)
for a_name,a in (("a=1.89 Gpc (3D mean spacing)",1.89),("a=2.03 Gpc (2D-matched)",2.03),("a=171 Mpc",0.171)):
    RES["_setup"][f"J_unit_{a_name}"]=G*MC**2/(a*GPC)
json.dump(RES,open("tri2d_results.json","w"),indent=1)
np.savez("tri2d_configs.npz",X=X,E=E,W=np.array(W),**{k:v for k,v in CONF.items()})
for k,o in RES.items():
    if not k.startswith("_"): print(f"{k:55s} like {o['like_links_per_node']:.3f} unlike {o['unlike_links_per_node']:.3f}  U_nn {o['U_nn_per_node']:+9.4f}  U_all {o['U_all_per_node']:+10.4f}")
# size scaling of the all-pairs energy (case A) for random / alternating rows / phase separation
SC={}
for lx,ly in ((22,24),(44,46),(66,70)):
    x,i_,j_,nb,e,w=lattice(lx,ly); p=phi_matrix(x,i_,j_,w,lx,ly); n=len(x)
    r=np.random.default_rng(1).permutation(n)<n//2
    for nm,c in (("random",r),("alternating rows",(j_%2==0)),("phase separation (best interface)",(i_<lx//2) if w[0]>=w[1] else (j_<ly//2))):
        m=np.where(c,1.0,-1.0); SC.setdefault(nm,[]).append((n,float(U_full(m,p)/n),float(U_nn(m,e)/n)))
    print(lx,ly,{k:v[-1] for k,v in SC.items()})
json.dump(SC,open("tri2d_size_scaling.json","w"),indent=1)
