"""Entropy / topology / toy energies of the foam s2 mesh faces, BEFORE vs AFTER, whole slab and bands by distance from the ring.
Same mesh and same kept faces as mesh_faces.py (functions reused verbatim). TOY ESTIMATE for the energies -- see notes in table.md."""
import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src=open("mesh_faces.py").read(); exec(src[:src.index("js=sys.argv[1]")])          # imports + faces(), valid(), AMAX
import seq_core
G=6.674e-11; C=2.998e8; V=3.0e5; GPC=seq_core.GPC
d=json.load(open("minedit_ent/b_foam_rim_s2_tidalcore_entropy.json")); base=d["tag"].replace("_entropy","")
Lb=EB.layout(base); La=EB.layout(base,d)
BANDS={"whole":(0,1e9),"0-3 Gpc":(0,3),"3-6 Gpc":(3,6),"6+ Gpc":(6,1e9)}
def dring2(p): return np.abs(np.hypot(p[...,0],p[...,1])-RMIN)
gam=None; out={}
for view,(X,r,v,M,_) in (("before",Lb),("after",La)):
    sl,E,T=MV.network(X,r,v); P=X[sl,:2]; X3=X[sl]*GPC; Ms=M[sl]; vs=v[sl]
    F,_=valid(P,faces(P,E)); A=np.array([a for _,a in F]); cen=np.array([P[f].mean(0) for f,_ in F]); dfc=dring2(cen)
    sides=np.array([len(f) for f,_ in F]); per=np.array([np.sum(np.linalg.norm(P[f]-np.roll(P[f],-1,0),axis=1)) for f,_ in F])*GPC
    Ll=np.linalg.norm(X3[E[:,0]]-X3[E[:,1]],axis=1); dlk=dring2(0.5*(P[E[:,0]]+P[E[:,1]]))
    Ulink=-G*Ms[E[:,0]]*Ms[E[:,1]]/Ll
    if gam is None:     # filament mass per length, fixed from BEFORE: half of each end cluster's mass per cluster-cluster link, over cc link length
        cc=~vs[E[:,0]]&~vs[E[:,1]]; deg=np.bincount(E[cc].ravel(),minlength=len(sl)); hm=np.where(deg>0,np.abs(Ms)/np.maximum(deg,1),0)
        gam=float(np.sum(hm[E[cc,0]]+hm[E[cc,1]])/np.sum(Ll[cc]))
        kT=float(np.mean(0.5*np.abs(Ms[~vs])*V**2))     # T_eff: mean kinetic energy per cluster body (k_B T_eff in J)
    R={}
    for b,(lo,hi) in BANDS.items():
        mf=(dfc>=lo)&(dfc<hi); ml=(dlk>=lo)&(dlk<hi); a=A[mf]; p=a/a.sum(); S=float(-np.sum(p*np.log(p))); N=int(mf.sum())
        sc=np.bincount(sides[mf]); ps=sc[sc>0]/sc.sum()
        R[b]=dict(n_faces=N,total_area_Gpc2=float(a.sum()),S_area=S,S_per_face_evenness=S/np.log(N),sides_mean=float(sides[mf].mean()),sides_mu2=float(sides[mf].var()),
                  S_sides=float(-np.sum(ps*np.log(ps))),sides_hist={int(k):int(c) for k,c in enumerate(sc) if c},n_links=int(ml.sum()),U_J=float(Ulink[ml].sum()),
                  perimeter_m=float(per[mf].sum()),Ewall_c2_J=float(gam*C**2*per[mf].sum()),Ewall_v2_J=float(gam*V**2*per[mf].sum()))
    out[view]=R
ch={}
for b in BANDS:
    B,A_=out["before"][b],out["after"][b]; c={k:A_[k]-B[k] for k in B if isinstance(B[k],(int,float))}
    c["U_frac"]=c["U_J"]/abs(B["U_J"]); c["Ewall_frac"]=c["Ewall_v2_J"]/B["Ewall_v2_J"]
    c["TdS_area_J"]=kT*c["S_area"]; c["TdS_evenness_J"]=kT*c["S_per_face_evenness"]
    c["dF_J"]=c["U_J"]+c["Ewall_v2_J"]-kT*c["S_area"]; ch[b]=c
res=dict(constants=dict(G=G,c=C,v=V,GPC_m=GPC,gamma_kg_per_m=gam,kT_eff_J=kT),before=out["before"],after=out["after"],change=ch,
  notes=["signed masses as in the model: clusters +M (kg), voids = negative mass deficit (DELTA_V*rho*V); so cluster-void links give U>0, void-void and cluster-cluster U<0",
         "U over mesh links only (not all pairs), 3D separations of the slab bodies; links assigned to a band by midpoint (projected) distance from the ring",
         "gamma fixed from BEFORE: each cluster's mass shared equally over its cluster-cluster links, summed, divided by total cc link length; E_wall = gamma*v^2*perimeter (v=300 km/s) and gamma*c^2*perimeter (rest-mass version)",
         "perimeter = sum of kept face perimeters in the slab projection (shared edges counted twice, as faces)",
         "T_eff: k_B*T_eff = mean (1/2) M v^2 over slab clusters, v=300 km/s; S is dimensionless (nats); dF = dU + dE_wall(v^2) - kT_eff*dS_area. TOY ESTIMATE."])
json.dump(res,open("mesh_energy.json","w"),indent=1)
# table
f_=lambda x: f"{x:.3e}" if abs(x)>=1e4 or (abs(x)<1e-2 and x!=0) else f"{x:.3f}"
rows=[("faces","n_faces"),("total area (Gpc²)","total_area_Gpc2"),("area entropy S (nats)","S_area"),("evenness S/ln N","S_per_face_evenness"),("mean sides","sides_mean"),
      ("sides variance μ2","sides_mu2"),("side-count entropy","S_sides"),("links","n_links"),("U over links (J)","U_J"),("perimeter (m)","perimeter_m"),
      ("E_wall γv²L (J)","Ewall_v2_J"),("E_wall γc²L (J)","Ewall_c2_J")]
md=["# foam s2 mesh: entropy, topology and toy energies (before → after)","",f"γ = {gam:.3e} kg/m (fixed from BEFORE); k_B·T_eff = {kT:.3e} J per body (½Mv², v = 300 km/s); G = 6.674e-11.",""]
for b in BANDS:
    md+= [f"## {b}","","| quantity | before | after | change |","|---|---|---|---|"]
    for lab,k in rows: md.append(f"| {lab} | {f_(out['before'][b][k])} | {f_(out['after'][b][k])} | {f_(ch[b][k])} |")
    c=ch[b]; md+= [f"| ΔU / |U| | | | {100*c['U_frac']:+.1f}% |",f"| ΔE_wall / E_wall | | | {100*c['Ewall_frac']:+.1f}% |",f"| T_eff·ΔS (area) (J) | | | {f_(c['TdS_area_J'])} |",
                   f"| ΔF = ΔU + ΔE_wall(v²) − T_eff·ΔS (J), toy | | | {f_(c['dF_J'])} |",""]
md+= ["## Notes",""]+["- "+n for n in res["notes"]]
open("table.md" if False else "mesh_energy_table.md","w").write("\n".join(md)+"\n")
# chart: 4 small panels, change by band
fig,ax=plt.subplots(1,5,figsize=(18,4.8)); fig.subplots_adjust(wspace=0.45,bottom=0.28,top=0.76)
bl=list(BANDS); cols=["#34495e","#8e44ad","#e67e22","#16a085"]; x=np.arange(4)
items=[("ΔS area entropy (nats)",[ch[b]["S_area"] for b in bl]),("ΔS evenness (S/ln N)",[ch[b]["S_per_face_evenness"] for b in bl]),
       ("ΔU over links (% of |U| before)",[100*ch[b]["U_frac"] for b in bl]),("ΔE_wall (% of before)",[100*ch[b]["Ewall_frac"] for b in bl]),
       ("ΔF toy (% of |U| before)",[100*ch[b]["dF_J"]/abs(out["before"][b]["U_J"]) for b in bl])]
for a_,(t,vals) in zip(ax,items):
    a_.bar(x,vals,color=cols); a_.axhline(0,color="0.3",lw=0.8); a_.set_xticks([]); a_.set_title(t.replace(" (","\n("),fontsize=12); a_.tick_params(labelsize=11)
    for xi,vv in zip(x,vals): a_.text(xi,vv,f"{vv:+.3g}",ha="center",va="bottom" if vv>=0 else "top",fontsize=10)
    lo,hi=min(0,min(vals)),max(0,max(vals)); pad=0.18*(hi-lo or 1); a_.set_ylim(lo-pad,hi+pad)
fig.suptitle("foam s2 mesh, AFTER − BEFORE by distance from the ring (toy energies)",fontsize=13,fontweight="bold",y=0.99)
sgn="lower" if ch["whole"]["dF_J"]<0 else "higher"
fig.text(0.5,0.885,f"ΔF (toy, whole slab) = {ch['whole']['dF_J']:.2e} J  →  free energy {sgn} after the edits;  T_eff·ΔS = {ch['whole']['TdS_area_J']:.2e} J",ha="center",fontsize=11)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=c,label=b) for c,b in zip(cols,bl)],loc="upper center",bbox_to_anchor=(0.5,0.14),ncol=4,frameon=False,fontsize=12)
fig.savefig("mesh_foam_s2_energy.png",dpi=140,bbox_inches="tight")
print(json.dumps({b:{k:(round(v_,4) if abs(v_)<1e4 else float("%.3e"%v_)) for k,v_ in ch[b].items()} for b in bl},indent=0)); print(gam,kT)
for b in bl: print(b,{k:out[v_][b][k] for v_ in ("before","after") for k in ("S_area","S_per_face_evenness","sides_mean","sides_mu2","S_sides","U_J","Ewall_v2_J")})
