import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
R=json.load(open("tri2d_results.json")); Z=np.load("tri2d_configs.npz"); X,E,W=Z["X"],Z["E"],Z["W"]; FS=12; FT=13; FL=11
# ---- configurations
def draw(ax,c,title):
    D=X[E[:,1]]-X[E[:,0]]; ok=(np.abs(D[:,0])<2)&(np.abs(D[:,1])<2)   # skip periodic wrap links
    ax.add_collection(LineCollection(np.stack([X[E[ok,0]],X[E[ok,1]]],1),colors="0.75",lw=0.35,zorder=1))
    ax.scatter(X[c,0],X[c,1],s=7,c="#c0392b",lw=0,zorder=2); ax.scatter(X[~c,0],X[~c,1],s=7,c="#2e86c1",lw=0,zorder=2)
    ax.set_aspect("equal"); ax.set_xlim(-0.8,W[0]); ax.set_ylim(-0.8,W[1]); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title,fontsize=FT)
gs=Z["MC anneal, all pairs (Ewald) [A]"]
fig,ax=plt.subplots(3,1,figsize=(8.2,22)); fig.subplots_adjust(hspace=0.22,bottom=0.06,top=0.97)
k=lambda n: R["A | "+n]
for a_,c,n,t in ((ax[0],gs,"MC anneal, all pairs (Ewald) [A]","(a) ground state: full phase separation"),(ax[1],Z["random 50/50"],"random 50/50","(b) random 50/50"),
                 (ax[2],Z["alternating rows (= best alternation)"],"alternating rows (= best alternation)","(c) alternating rows (best alternation)")):
    o=k(n); draw(a_,c,f"{t}\nlike / unlike neighbours per node: {o['like_links_per_node']:.2f} / {o['unlike_links_per_node']:.2f}\nenergy per node: {o['U_nn_per_node']:+.2f} (neighbours), {o['U_all_per_node']:+.2f} (all pairs)")
fig.legend(handles=[Line2D([],[],marker="o",ls="",ms=6,color="#c0392b",label="cluster (+M)"),Line2D([],[],marker="o",ls="",ms=6,color="#2e86c1",label="void (−M)"),
    Line2D([],[],color="0.6",lw=1,label="link (6 per node, length a)")],loc="lower center",ncol=3,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.005),
    title="2,024-site triangular lattice, periodic box, equal numbers\nenergy in units of G·M²/a, case A (|M_void| = M_cluster)",title_fontsize=FL)
fig.savefig("tri2d_configs.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# ---- bar chart (energy per node relative to random)
names=["alternating rows (= best alternation)","double-row stripes","small clumps","big clumps","phase separation, row interface","phase separation, column interface",
       "MC anneal, nearest-neighbour links","MC anneal, all pairs (Ewald)"]
lab=["alternating rows\n(best alternation)","double-row stripes","small clumps","big clumps","phase separation,\nrow interface","phase separation,\ncolumn interface","MC anneal\n(neighbour links)","MC anneal\n(all pairs)"]
fig,ax=plt.subplots(2,1,figsize=(10,12)); fig.subplots_adjust(hspace=0.32,bottom=0.1)
for a_,case,tt in ((ax[0],"A","(A) equal magnitudes, |M_void| = M_cluster"),(ax[1],"B","(B) web values, M_cluster = 1.40e50 kg, M_void = −8.39e50 kg")):
    g=lambda n: R[f"{case} | {n}"+(f" [{case}]" if n.startswith("MC") else "")]
    r0=R[f"{case} | random 50/50"]; y=np.arange(len(names))
    nn=[g(n)["U_nn_per_node"]-r0["U_nn_per_node"] for n in names]; ap=[g(n)["U_all_per_node"]-r0["U_all_per_node"] for n in names]
    a_.barh(y-0.2,nn,0.38,color="#7f8c8d",label="nearest-neighbour links only"); a_.barh(y+0.2,ap,0.38,color="#8e44ad",label="all pairs (Ewald)")
    a_.axvline(0,color="k",lw=0.8); a_.set_yticks(y); a_.set_yticklabels(lab,fontsize=10.5); a_.invert_yaxis(); a_.grid(axis="x",alpha=0.3)
    a_.set_title(tt,fontsize=FT,loc="left"); a_.set_xlabel("energy per node minus random, in G·M²_cluster / a   (← more bound)",fontsize=FS); a_.tick_params(axis="x",labelsize=10.5)
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",ncol=2,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.02))
fig.suptitle("Energy per node by configuration (2,024-site triangular lattice)",fontsize=FT+1,fontweight="bold",y=0.97)
fig.savefig("tri2d_energy_bars.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# ---- size scaling
S=json.load(open("tri2d_size_scaling.json")); fig,ax=plt.subplots(figsize=(8.5,6.2)); fig.subplots_adjust(bottom=0.25)
for nm,c in (("phase separation (best interface)","#c0392b"),("random","#e67e22"),("alternating rows","#2e86c1")):
    v=np.array(S[nm]); ax.plot(np.sqrt(v[:,0]),v[:,1],"o-",color=c,label=nm+", all pairs"); ax.plot(np.sqrt(v[:,0]),v[:,2],"s--",color=c,alpha=0.6,label=nm+", neighbours only")
ax.set_xlabel("box size, √(number of sites)",fontsize=FS); ax.set_ylabel("energy per node (G·M²/a), case A",fontsize=FS); ax.grid(alpha=0.3); ax.tick_params(labelsize=10.5)
ax.set_title("All-pairs energy of the separated state keeps falling as the box grows",fontsize=FT)
fig.legend(loc="lower center",ncol=2,frameon=False,fontsize=10,bbox_to_anchor=(0.5,0.0)); fig.savefig("tri2d_size_scaling.png",dpi=140,bbox_inches="tight"); plt.close(fig)
print("ok")
