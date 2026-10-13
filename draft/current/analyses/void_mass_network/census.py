"""Full link census of one generated network, before and after the non-crossing rule (planar.py):
cluster-cluster (filaments: direct foam edges 'cc', and straight cluster-cluster links through one empty web node 'cc2'),
cluster-void (bipartite 'cv'), void-void (shared walls 'vv').
Census region: nodes / link midpoints inside R_count of a sub-ball R_sub where the crossing search is done.
Writes census_<Mmin>.json, results/census_arrays.npz."""
import os,sys,json,numpy as np,scipy.sparse as sp
ME=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,ME); import vmn_gen as G, planar as PL
Mmin=float(sys.argv[1]) if len(sys.argv)>1 else 3e15; RT=float(sys.argv[2]) if len(sys.argv)>2 else 0.05
R_SUB,R_CNT=7000.0,4500.0
g=G.generate(Mmin=Mmin,seed=7); S=g["S"]; XC=g["XC"]; MC=g["MC"]; Reff=g["Reff"]; ncl=len(MC); nv=len(S)
isC=g["node_cl"]>=0; fe=g["foam_edges"]
cc=np.sort(g["node_cl"][fe[isC[fe[:,0]]&isC[fe[:,1]]]],1)
nn=len(g["npos"]); Adj=sp.coo_matrix((np.ones(2*len(fe)),(np.r_[fe[:,0],fe[:,1]],np.r_[fe[:,1],fe[:,0]])),shape=(nn,nn)).tocsr()
Cn=np.where(isC)[0]; En=np.where(~isC)[0]; ACE=Adj[Cn][:,En]; two=(ACE@ACE.T).tocoo(); m=two.row<two.col
cc2=np.sort(np.c_[g["node_cl"][Cn[two.row[m]]],g["node_cl"][Cn[two.col[m]]]],1); cc2=np.unique(cc2,axis=0)
sA={tuple(x) for x in cc}; cc2=np.array([x for x in map(tuple,cc2) if x not in sA])
r=PL.planarize(g,cc,cc2,R_sub=R_SUB,rt_frac=RT); E,T,L,kept,why,pos=r["E"],r["T"],r["L"],r["kept"],r["why"],r["pos"]
NAMES=["cc","cc2","cv","vv"]
Rn=np.r_[np.zeros(ncl),Reff]; isvoid=np.r_[np.zeros(ncl,bool),np.ones(nv,bool)]
nodein=np.linalg.norm(pos,axis=1)<R_CNT; mid=0.5*(pos[E[:,0]]+pos[E[:,1]]); midin=np.linalg.norm(mid,axis=1)<R_CNT
# Rv attached to each link: cv -> its void; vv -> mean of the two; cc/cc2 -> mean Reff of voids the two clusters share
cvl=g["links"]; Binc=sp.csr_matrix((np.ones(len(cvl)),(cvl[:,0],cvl[:,1])),shape=(ncl,nv)); BR=Binc.multiply(Reff[None,:]).tocsr()
Rl=np.full(len(E),np.nan)
for t in (0,1):
    mk=T==t; a,b=E[mk,0],E[mk,1]
    if mk.sum()==0: continue
    num=np.asarray(BR[a].multiply(Binc[b]).sum(1)).ravel(); den=np.asarray(Binc[a].multiply(Binc[b]).sum(1)).ravel()
    Rl[mk]=np.where(den>0,num/np.maximum(den,1),np.nan)
    # cc2 pairs may share no void: use mean Reff of all voids of both clusters
    nb=den==0
    if nb.any():
        ra=np.asarray(BR[a[nb]].sum(1)).ravel()/np.asarray(Binc[a[nb]].sum(1)).ravel(); rb=np.asarray(BR[b[nb]].sum(1)).ravel()/np.asarray(Binc[b[nb]].sum(1)).ravel()
        tmp=Rl[mk]; tmp[nb]=0.5*(ra+rb); Rl[mk]=tmp
mk=T==2; Rl[mk]=Rn[E[mk,1]]; mk=T==3; Rl[mk]=0.5*(Rn[E[mk,0]]+Rn[E[mk,1]])
mt=np.quantile(MC,[1/3,2/3]); mbin=np.digitize(MC,mt); vt=np.quantile(Reff[g["void_in"]],[1/3,2/3]); vbin=np.digitize(Reff,vt)
def st(x):
    x=np.asarray(x,float); x=x[np.isfinite(x)]
    if len(x)==0: return dict(n=0)
    return dict(mean=float(x.mean()),sd=float(x.std()),median=float(np.median(x)),p16=float(np.percentile(x,16)),p84=float(np.percentile(x,84)),min=float(x.min()),max=float(x.max()),n=int(len(x)))
def hist(x): b=np.bincount(np.asarray(x,int)); return {int(i):int(c) for i,c in enumerate(b) if c}
def angles(sel_nodes,Em,maxn=6000,seed=0):
    rng=np.random.default_rng(seed); E2=np.vstack([Em,Em[:,::-1]]); o=np.argsort(E2[:,0]); E2=E2[o]; stt=np.searchsorted(E2[:,0],np.arange(len(pos)+1))
    nodes=sel_nodes if len(sel_nodes)<=maxn else rng.choice(sel_nodes,maxn,replace=False); out=[]
    for c in nodes:
        P=pos[E2[stt[c]:stt[c+1],1]]
        if len(P)<2: continue
        u=P-pos[c]; u/=np.linalg.norm(u,axis=1)[:,None]; d=np.clip(u@u.T,-1,1); i,j=np.triu_indices(len(P),1); out.append(np.degrees(np.arccos(d[i,j])))
    return np.concatenate(out) if out else np.array([])
cl_in=np.where(nodein[:ncl])[0]; v_in=ncl+np.where(nodein[ncl:])[0]
res={}; arrays={}
for setname,sel in (("before",np.ones(len(E),bool)),("after",kept)):
    R={}
    deg={t:np.bincount(E[sel&(T==t)].ravel(),minlength=len(pos)) for t in range(4)}
    ccall=deg[0]+deg[1]
    R["filaments_per_cluster_direct"]=dict(stats=st(deg[0][cl_in]),hist=hist(deg[0][cl_in]),frac_zero=float((deg[0][cl_in]==0).mean()),by_mass={f"third_{k}":st(deg[0][cl_in][mbin[cl_in]==k]) for k in range(3)})
    R["filaments_per_cluster_all"]=dict(stats=st(ccall[cl_in]),hist=hist(ccall[cl_in]),frac_zero=float((ccall[cl_in]==0).mean()),by_mass={f"third_{k}":st(ccall[cl_in][mbin[cl_in]==k]) for k in range(3)})
    R["voids_per_cluster"]=dict(stats=st(deg[2][cl_in]),hist=hist(deg[2][cl_in]),by_mass={f"third_{k}":st(deg[2][cl_in][mbin[cl_in]==k]) for k in range(3)})
    vb=vbin[v_in-ncl]
    R["clusters_per_void"]=dict(stats=st(deg[2][v_in]),hist=hist(deg[2][v_in]),by_void_size={f"third_{k}":st(deg[2][v_in][vb==k]) for k in range(3)})
    R["neighbour_voids_per_void"]=dict(stats=st(deg[3][v_in]),hist=hist(deg[3][v_in]),by_void_size={f"third_{k}":st(deg[3][v_in][vb==k]) for k in range(3)})
    for t in range(4):
        mm=sel&(T==t)&midin; R[f"length_{NAMES[t]}_Mpc"]=st(L[mm]); R[f"length_{NAMES[t]}_over_Rv"]=st(L[mm]/Rl[mm])
        arrays[f"{setname}_L_{NAMES[t]}"]=L[mm]; arrays[f"{setname}_LR_{NAMES[t]}"]=(L/Rl)[mm]
    a_cc=angles(cl_in,E[sel&(T<=1)]); a_cv_c=angles(cl_in,E[sel&(T==2)]); a_vv=angles(v_in,E[sel&(T==3)]); a_cv_v=angles(v_in,E[sel&(T==2)])
    R["angles_deg"]=dict(filaments_at_cluster=st(a_cc),void_links_at_cluster=st(a_cv_c),wall_links_at_void=st(a_vv),cluster_links_at_void=st(a_cv_v))
    for k,v in (("deg_cc",deg[0][cl_in]),("deg_ccall",ccall[cl_in]),("deg_cv_c",deg[2][cl_in]),("deg_cv_v",deg[2][v_in]),("deg_vv",deg[3][v_in]),
                ("ang_cc",a_cc),("ang_cvc",a_cv_c),("ang_vv",a_vv),("ang_cvv",a_cv_v)): arrays[f"{setname}_{k}"]=v
    res[setname]=R
rem={}
for t in range(4):
    mm=(T==t)&midin; rem[NAMES[t]]=dict(candidates=int(mm.sum()),kept=int((kept&mm).sum()),removed_crossing=int((mm&(why==2)).sum()),removed_through_body=int((mm&(why==1)).sum()),
        removed_frac=float(1-(kept&mm).sum()/max(mm.sum(),1)))
arrays.update(mbin=mbin[cl_in],vbin=vbin[v_in-ncl],alldeg=np.bincount(fe.ravel(),minlength=nn)[g["node_in"]])
out=dict(Mmin=Mmin,tube_radius_Mpc=r["rt"],tube_radius_frac_of_median_Rv=RT,R_sub_Mpc=R_SUB,R_count_Mpc=R_CNT,n_clusters_counted=int(len(cl_in)),n_voids_counted=int(len(v_in)),
         cluster_fraction_of_web_nodes=float(ncl/g["n_nodes_in"]),mass_third_edges_Msun=mt.tolist(),void_Reff_third_edges_Mpc=vt.tolist(),
         all_web_nodes_degree=st(arrays["alldeg"]),removed=rem,before=res["before"],after=res["after"])
tag=f"{Mmin:.0e}" if RT==0.05 else f"{Mmin:.0e}_rt{RT}"
json.dump(out,open(os.path.join(ME,f"census_{tag}.json"),"w"),indent=1)
if RT==0.05:
    np.savez_compressed(os.path.join(ME,"results","census_arrays.npz"),**arrays)
    np.savez_compressed(os.path.join(ME,"results","planar_links.npz"),E=E,T=T,kept=kept,why=why,pos=pos,Rv_core=np.r_[np.zeros(ncl),g["Rv"]],MC=MC)
print(json.dumps(rem,indent=1))
for s in ("before","after"):
    R=res[s]; print(s,{k:round(v["stats"]["mean"],2) for k,v in R.items() if isinstance(v,dict) and "stats" in v})
    print("  len med",{k:round(v["median"],1) for k,v in R.items() if k.startswith("length") and "median" in v})
    print("  ang med",{k:round(v["median"],1) for k,v in R["angles_deg"].items() if "median" in v})
