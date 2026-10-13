import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
C=json.load(open("counts.json")); S=[s for s in C["scenarios"] if s.get("reach",True) and "M_plus" in s]
def best(geom,sym):
    c=[s for s in S if (s["geom"]==geom or geom=="both") and s["sym"]==sym and s["left"]<=1.001 and (s.get("dmin") or 0)>=-1]
    if not c: c=[s for s in S if s["geom"]==geom and s["sym"]==sym and s["left"]<=1.001]
    return min(c,key=lambda s:s["M_plus"]) if c else None
def geom_note(b,g): return False
cats=[("inside","mirror"),("inside","no-mirror"),("outside","mirror"),("outside","no-mirror"),("both","mirror"),("both","no-mirror")]
fig,ax=plt.subplots(figsize=(13,7.5)); x=np.arange(len(cats)); w=0.36
for i,(g,sy) in enumerate(cats):
    b=best(g,sy)
    if b is None: ax.text(i,2e5,"no 1× fit",ha="center",fontsize=11,rotation=90); continue
    if geom_note(b,g): pass
    ax.bar(i-w/2,b["n_clusters"],w,color="#c0392b"); ax.bar(i+w/2,b["n_voids"],w,color="#2471a3")
    NATg=C["natural"][g]; cl=C["constants"]
    for dx,nat in ((-w/2,NATg["clu_mass"]/cl["rich_cluster_kg"]),(w/2,NATg["void_def"]/cl["void_deficit_kg"])): ax.plot([i+dx-w/2,i+dx+w/2],[nat,nat],"k-",lw=3)
    flag="" if (b.get("dmin") or 0)>=-1 else "\nδ<−1: unphysical"
    ax.text(i,max(b["n_clusters"],b["n_voids"])*1.6,("field" if b["kind"]=="field" else "points")+(" (all inside)" if g=="both" and b["geom"]=="inside" else "")+flag,ha="center",fontsize=11,color="0.3")
ax.set_yscale("log"); ax.set_ylim(1e5,3e10)
ax.set_xticks(x); ax.set_xticklabels([f"{g}\n{'mirror allowed' if s=='mirror' else 'no mirror'}" for g,s in cats],fontsize=12)
ax.set_ylabel("number of objects needed (extra, beyond normal)")
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#c0392b",label="rich clusters (1e15 M☉) for the dense side"),Patch(color="#2471a3",label="voids (R 29 Mpc, δ −0.4) for the empty side"),
                   Line2D([0],[0],color="k",lw=3,label="how many that region holds anyway")],loc="upper left",frameon=False)
ax.set_title("Realistic objects needed to bring the CMB leftover to 1× (best fit per case)")
fig.tight_layout(); fig.savefig("objects_needed_by_scenario.png",dpi=110); plt.close(fig)

# map: inside, mirror-allowed smooth field with per-lobe counts
L=json.load(open("lobes_inside.json")); z=np.load("field_inside_mirror.npz"); d=z["arr_0"]; P=z["arr_1"]; H=0.3
sel=np.abs(P[:,2])<H/2+1e-6
fig,ax=plt.subplots(figsize=(10.5,10))
sc=ax.scatter(P[sel,0],P[sel,1],c=d[sel],cmap="RdBu_r",vmin=-0.05,vmax=0.05,s=9,marker="s",linewidths=0)
ax.add_patch(Circle((0,0),14.3,fill=False,ls="--",color="k",lw=1.5)); ax.plot(0,0,"k+",ms=16,mew=2.5)
for l in L:
    a=np.radians(l["lon"]); rr=7.2
    txt=(f"+{l['n_clusters']:.1e} clusters\n(or {l['n_superclusters']:.0e} superclusters)" if l["sign"]=="+" else f"{l['n_voids']:.1e} voids\n(or {l['n_supervoids']:.0e} supervoids)")
    txt=txt.replace("e+0","e").replace("e+","e")
    ax.text(rr*np.cos(a),rr*np.sin(a),txt,ha="center",va="center",fontsize=10.5,color="#922b21" if l["sign"]=="+" else "#1a5276",
            bbox=dict(facecolor="white",alpha=0.85,edgecolor="none",pad=1.5))
ax.set_xlim(-16,16); ax.set_ylim(-16,16); ax.set_aspect("equal"); ax.set_xlabel("Gpc"); ax.set_ylabel("Gpc")
cb=fig.colorbar(sc,ax=ax,shrink=0.7,pad=0.02); cb.set_label("δ  (red over-dense, blue under-dense)")
ax.set_title("Best case: inside our sphere, smooth lobes near the CMB ring\n(each lobe ≈ 400–700 Gpc³, 10–11 Gpc from us; counts = extra objects per lobe)",fontsize=13)
fig.tight_layout(); fig.savefig("map_where_objects_go.png",dpi=110); plt.close(fig)
