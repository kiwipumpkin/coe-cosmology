"""Entropy / topology / toy energies for BOTH meshes of foam s2 (Chris, 10 Oct):
  'faces' = non-crossing network faces (mesh_view.network + mesh_faces.faces/valid, same kept faces as mesh_foam_s2_faces.json)
  'tri'   = Delaunay triangulation of the same slab nodes (tri_view.tri, triangles with an edge > 6 Gpc dropped)
Bands by distance of the face centroid (links: midpoint) from the ring: whole, 0-3, 3-6, 6+ Gpc.
Gap handling: 'all' = every kept face/link in the band; 'common' = only faces/links whose centroid/midpoint lies in area covered by a face
of that mesh in BOTH views (so newly covered gap area next to the ring is not counted). TOY ESTIMATE for the energies."""
import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.path import Path
src=open("mesh_faces.py").read(); exec(src[:src.index("js=sys.argv[1]")])
src=open("tri_view.py").read(); exec(src[src.index("LMAX="):src.index("def links(")])
import seq_core
from scipy.spatial import Delaunay
G=6.674e-11; C=2.998e8; V=3.0e5; GPC=seq_core.GPC
d=json.load(open("minedit_ent/b_foam_rim_s2_tidalcore_entropy.json")); base=d["tag"].replace("_entropy","")
LAY={"before":EB.layout(base),"after":EB.layout(base,d)}
BANDS={"whole":(0,1e9),"0-3 Gpc":(0,3),"3-6 Gpc":(3,6),"6+ Gpc":(6,1e9)}
def dring2(p): return np.abs(np.hypot(p[...,0],p[...,1])-RMIN)
def build(kind,view):
    X,r,v,M,_=LAY[view]
    if kind=="faces":
        sl,E,_=MV.network(X,r,v); P=X[sl,:2]; F=[np.array(f) for f,_ in valid(P,faces(P,E))[0]]
    else:
        sl,P,S,E,_=tri(X,r); F=[s for s in S]
    return dict(sl=sl,P=P,E=E,F=F,X3=X[sl]*GPC,M=M[sl],v=v[sl])
g=np.linspace(-30,30,601)
def cover(m):
    out=np.zeros((len(g),len(g)),bool); GX,GY=np.meshgrid(g,g); pts=np.c_[GX.ravel(),GY.ravel()]; o=out.ravel()
    for f in m["F"]:
        q=m["P"][f]; lo,hi=q.min(0),q.max(0); s=np.where(np.all((pts>=lo)&(pts<=hi),axis=1))[0]
        if len(s): o[s[Path(q).contains_points(pts[s])]]=True
    return o.reshape(out.shape)
def incov(cv,p): i=np.clip(np.round((p[:,1]+30)/0.1).astype(int),0,600); j=np.clip(np.round((p[:,0]+30)/0.1).astype(int),0,600); return cv[i,j]
kT=None; RES={}
for kind in ("faces","tri"):
    m={vw:build(kind,vw) for vw in ("before","after")}; cv={vw:cover(m[vw]) for vw in m}; both=cv["before"]&cv["after"]
    b=m["before"]; cc=~b["v"][b["E"][:,0]]&~b["v"][b["E"][:,1]]; E=b["E"]; Ll=np.linalg.norm(b["X3"][E[:,0]]-b["X3"][E[:,1]],axis=1)
    deg=np.bincount(E[cc].ravel(),minlength=len(b["sl"])); hm=np.where(deg>0,np.abs(b["M"])/np.maximum(deg,1),0)
    gam=float(np.sum(hm[E[cc,0]]+hm[E[cc,1]])/np.sum(Ll[cc]))
    if kT is None: kT=float(np.mean(0.5*np.abs(b["M"][~b["v"]])*V**2))
    out={}
    for vw,mm in m.items():
        P=mm["P"]; A=np.array([0.5*abs(np.sum(P[f,0]*np.roll(P[f,1],-1)-np.roll(P[f,0],-1)*P[f,1])) for f in mm["F"]])
        cen=np.array([P[f].mean(0) for f in mm["F"]]); sides=np.array([len(f) for f in mm["F"]])
        per=np.array([np.sum(np.linalg.norm(P[f]-np.roll(P[f],-1,0),axis=1)) for f in mm["F"]])*GPC
        E=mm["E"]; mid=0.5*(P[E[:,0]]+P[E[:,1]]); Ll=np.linalg.norm(mm["X3"][E[:,0]]-mm["X3"][E[:,1]],axis=1); U=-G*mm["M"][E[:,0]]*mm["M"][E[:,1]]/Ll
        fc,lc=incov(both,cen),incov(both,mid); dfc,dlk=dring2(cen),dring2(mid)
        for gp,(fm,lm) in (("all",(np.ones(len(A),bool),np.ones(len(U),bool))),("common",(fc,lc))):
            for bn,(lo,hi) in BANDS.items():
                mf=fm&(dfc>=lo)&(dfc<hi); ml=lm&(dlk>=lo)&(dlk<hi); a=A[mf]; p=a/a.sum(); S=float(-np.sum(p*np.log(p))); N=int(mf.sum())
                sc=np.bincount(sides[mf]); ps=sc[sc>0]/sc.sum()
                out.setdefault(gp,{}).setdefault(bn,{})[vw]=dict(n_faces=N,area_Gpc2=float(a.sum()),S_area=S,evenness=S/np.log(N),sides_mean=float(sides[mf].mean()),
                    sides_mu2=float(sides[mf].var()),S_sides=float(-np.sum(ps*np.log(ps))),n_links=int(ml.sum()),U_J=float(U[ml].sum()),perimeter_m=float(per[mf].sum()),
                    Ewall_v2_J=float(gam*V**2*per[mf].sum()),Ewall_c2_J=float(gam*C**2*per[mf].sum()))
    for gp in out:
        for bn in out[gp]:
            B,A_=out[gp][bn]["before"],out[gp][bn]["after"]; c={k:A_[k]-B[k] for k in B}
            c.update(U_frac=c["U_J"]/abs(B["U_J"]),Ewall_frac=c["Ewall_v2_J"]/B["Ewall_v2_J"],TdS_J=kT*c["S_area"],
                     dF_v2_J=c["U_J"]+c["Ewall_v2_J"]-kT*c["S_area"],dF_c2_J=c["U_J"]+c["Ewall_c2_J"]-kT*c["S_area"]); out[gp][bn]["change"]=c
    RES[kind]=dict(gamma_kg_per_m=gam,coverage_frac_before=float(cv["before"].mean()),coverage_frac_after=float(cv["after"].mean()),bands=out)
RES["constants"]=dict(G=G,c=C,v=V,kT_eff_J=kT,note="signed masses as in the model: clusters +M (kg); voids = negative mass deficit (DELTA_V*rho*V, about -8.4e50 kg); so cluster-void links give U>0. U over mesh links only, 3D separations. gamma per mesh from its BEFORE cluster-cluster links (each cluster's mass shared over its cc links), held fixed. E_wall = gamma*v^2*perimeter (v=300 km/s) and gamma*c^2*perimeter (upper bound). k_B T_eff = mean 1/2 M v^2 of slab clusters. dF = dU + dE_wall - kT_eff dS_area. TOY ESTIMATE.")
json.dump(RES,open("mesh_energy_compare.json","w"),indent=1)
# ---- table
f_=lambda x: f"{x:.2e}" if abs(x)>=1e4 or (abs(x)<1e-2 and x!=0) else f"{x:.3f}"
Q=[("faces","n_faces"),("area (Gpc²)","area_Gpc2"),("area entropy S","S_area"),("evenness S/lnN","evenness"),("mean sides","sides_mean"),("sides μ2","sides_mu2"),("side-count entropy","S_sides"),
   ("links","n_links"),("U (J)","U_J"),("E_wall γv²L (J)","Ewall_v2_J"),("E_wall γc²L (J)","Ewall_c2_J")]
md=["# foam s2: entropy, topology, toy energies — non-crossing face mesh vs Delaunay triangulation","",
    f"k_B·T_eff = {kT:.3e} J (½Mv², v = 300 km/s). γ: faces {RES['faces']['gamma_kg_per_m']:.3e} kg/m, tri {RES['tri']['gamma_kg_per_m']:.3e} kg/m.",
    "Sides for the triangulation are 3 by construction (μ2 = 0, side-count entropy = 0): trivial, shown for completeness.",
    "'all' = every kept face/link in the band; 'common' = only those whose centre lies in area covered in BOTH views (newly covered gap area excluded).",""]
rows=[("all",b) for b in BANDS]+[("common","0-3 Gpc"),("common","whole")]
for kind in ("faces","tri"):
    md+= [f"## {'non-crossing face mesh' if kind=='faces' else 'Delaunay triangulation'}","","| band | "+" | ".join(q for q,_ in Q)+" | ΔU/|U| | ΔE_wall/E | T_eff·ΔS (J) | ΔF v² (J) | ΔF c² (J) |","|---"*(len(Q)+6)+"|"]
    for gp,bn in rows:
        o=RES[kind]["bands"][gp][bn]; c=o["change"]
        md.append(f"| {bn}{' (common area)' if gp=='common' else ''} | "+" | ".join(f"{f_(o['before'][k])} → {f_(o['after'][k])}" for _,k in Q)+
                  f" | {100*c['U_frac']:+.1f}% | {100*c['Ewall_frac']:+.1f}% | {f_(c['TdS_J'])} | {f_(c['dF_v2_J'])} | {f_(c['dF_c2_J'])} |")
    md.append("")
md+= ["## Notes","","- "+RES["constants"]["note"]]
open("mesh_energy_compare_table.md","w").write("\n".join(md)+"\n")
# ---- paired bar chart: change by band, faces vs tri
labs=["whole","0-3 Gpc","0-3 Gpc\n(common area)","3-6 Gpc","6+ Gpc"]; keys=[("all","whole"),("all","0-3 Gpc"),("common","0-3 Gpc"),("all","3-6 Gpc"),("all","6+ Gpc")]
items=[("ΔS area entropy (nats)",lambda c,o:c["S_area"]),("Δ evenness S/ln N",lambda c,o:c["evenness"]),("ΔU over links\n(% of |U| before)",lambda c,o:100*c["U_frac"]),
       ("ΔE_wall\n(% of before)",lambda c,o:100*c["Ewall_frac"]),("ΔF toy, v² walls\n(% of |U| before)",lambda c,o:100*c["dF_v2_J"]/abs(o["before"]["U_J"]))]
fig,ax=plt.subplots(1,5,figsize=(19,5.2)); fig.subplots_adjust(wspace=0.38,bottom=0.3,top=0.82); x=np.arange(len(keys))
for a_,(t,fn) in zip(ax,items):
    for o_,kind,col in ((-0.2,"faces","#8e44ad"),(0.2,"tri","#16a085")):
        vals=[fn(RES[kind]["bands"][gp][bn]["change"],RES[kind]["bands"][gp][bn]) for gp,bn in keys]; a_.bar(x+o_,vals,0.38,color=col)
    a_.axhline(0,color="0.3",lw=0.8); a_.set_xticks(x); a_.set_xticklabels(labs,fontsize=9.5,rotation=40,ha="right"); a_.set_title(t,fontsize=12); a_.tick_params(axis="y",labelsize=10.5); a_.grid(axis="y",alpha=0.3)
from matplotlib.patches import Patch
fig.suptitle("foam s2, AFTER − BEFORE: non-crossing face mesh vs Delaunay triangulation (toy energies)",fontsize=13,fontweight="bold",y=0.97)
fig.legend(handles=[Patch(color="#8e44ad",label="non-crossing face mesh"),Patch(color="#16a085",label="Delaunay triangulation")],loc="upper center",bbox_to_anchor=(0.5,0.06),ncol=2,frameon=False,fontsize=12)
fig.savefig("mesh_foam_s2_energy_compare.png",dpi=140,bbox_inches="tight")
for kind in ("faces","tri"):
    print(kind,RES[kind]["gamma_kg_per_m"])
    for gp,bn in rows:
        o=RES[kind]["bands"][gp][bn]; c=o["change"]; B=o["before"]; A_=o["after"]
        print(f" {gp:6s} {bn:8s} N {B['n_faces']}->{A_['n_faces']} S {B['S_area']:.3f}->{A_['S_area']:.3f} ev {B['evenness']:.3f}->{A_['evenness']:.3f} mu2 {B['sides_mu2']:.2f}->{A_['sides_mu2']:.2f} Ss {B['S_sides']:.3f}->{A_['S_sides']:.3f} dU {c['U_J']:.2e} ({100*c['U_frac']:+.1f}%) dEw {100*c['Ewall_frac']:+.1f}% v2 {c['Ewall_v2_J']:.2e} c2 {c['Ewall_c2_J']:.2e} TdS {c['TdS_J']:.2e} dFv2 {c['dF_v2_J']:.2e} dFc2 {c['dF_c2_J']:.2e}")
