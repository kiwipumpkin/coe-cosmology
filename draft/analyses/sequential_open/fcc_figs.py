import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
D="fcc_lattice/"; R=json.load(open(D+"ladder.json")); RW=np.load(D+"ladder_raw.npz")
FS=12; FT=13; FL=11
raw=lambda k,q: RW[f"{k}|{q}"]
# ---------- 1) unit cell
a=171.0; fig=plt.figure(figsize=(8,8.4)); ax=fig.add_subplot(111,projection="3d"); ax.set_proj_type("ortho")
for s in [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]: pass
c=np.array([[i,j,k] for i in (0,1) for j in (0,1) for k in (0,1)],float)
for p in c:
    for q in c:
        if np.sum(np.abs(p-q))==1 and tuple(p)<tuple(q): ax.plot(*np.c_[p,q]*a,color="0.6",lw=1)
fc=np.array([[.5,.5,0],[.5,.5,1],[.5,0,.5],[.5,1,.5],[0,.5,.5],[1,.5,.5]])
S=np.vstack([c,fc])*a; ax.scatter(*S.T,s=260,c="#c0392b",edgecolor="#7b241c",depthshade=False,zorder=5)
tv=np.array([.25,.25,.25])*a; tc=np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]])*a
for i in range(4):
    for j in range(i+1,4): ax.plot(*np.c_[tc[i],tc[j]],color="#2e86c1",lw=2.4)
    ax.plot(*np.c_[tv,tc[i]],color="#2e86c1",lw=0.8,ls=":")
ax.add_collection3d(Poly3DCollection([tc[[0,1,2]],tc[[0,1,3]],tc[[0,2,3]],tc[[1,2,3]]],fc="#aed6f1",alpha=0.25,ec="none"))
ax.scatter(*tv,s=330,c="#2e86c1",edgecolor="#1b4f72",depthshade=False)
ov=np.array([.5,.5,.5])*a; oc=fc*a
for i in range(6):
    for j in range(i+1,6):
        if abs(np.linalg.norm(oc[i]-oc[j])-a/np.sqrt(2))<1: ax.plot(*np.c_[oc[i],oc[j]],color="#8e44ad",lw=1.4,ls="--")
ax.scatter(*ov,s=330,facecolors="white",edgecolors="#8e44ad",linewidths=2.2,depthshade=False)
ax.set_xlabel("x (Mpc)",fontsize=FL); ax.set_ylabel("y (Mpc)",fontsize=FL); ax.set_zlabel("z (Mpc)",fontsize=FL)
for t in (ax.set_xticks,ax.set_yticks,ax.set_zticks): t([0,85.5,171])
ax.tick_params(labelsize=10); ax.view_init(elev=18,azim=-58); ax.set_box_aspect((1,1,1))
ax.set_title("One FCC unit cell, edge 171 Mpc (Einasto period)",fontsize=FT,fontweight="bold")
fig.legend(handles=[Line2D([],[],marker="o",ls="",ms=12,mfc="#c0392b",mec="#7b241c",label="cluster (FCC site: 8 corners + 6 face centres)"),
    Line2D([],[],marker="o",ls="",ms=13,mfc="#2e86c1",mec="#1b4f72",label="void in a tetrahedral hole"),
    Line2D([],[],color="#2e86c1",lw=2.4,label="its 4-cluster pyramid (tetrahedron)"),
    Line2D([],[],marker="o",ls="",ms=13,mfc="white",mec="#8e44ad",mew=2,label="octahedral hole (empty here), for comparison"),
    Line2D([],[],color="#8e44ad",ls="--",lw=1.4,label="its 6-cluster octahedron")],
    loc="lower center",ncol=1,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.0),title="8 tetrahedral holes and 4 octahedral holes per cell; one of each shown",title_fontsize=FL)
fig.subplots_adjust(bottom=0.24,top=0.95); fig.savefig(D+"fcc_cell.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# ---------- histograms
def hist3(series,fn,title):
    fig,ax=plt.subplots(3,1,figsize=(9,12)); fig.subplots_adjust(hspace=0.42,bottom=0.14,top=0.93)
    vb=np.logspace(-3,1.3,70); lb=np.linspace(0,4,61); db=np.arange(3.5,40.5,1)
    for k,lab,col,lw,fill in series:
        for a_,q,b in ((ax[0],"vol",vb),(ax[2],"L",lb),(ax[1],"deg",db)):
            x=raw(k,q); w=np.ones(len(x))/len(x)
            if fill: a_.hist(x,b,weights=w,color=col,alpha=0.35,label=lab)
            else: a_.hist(x,b,weights=w,histtype="step",color=col,lw=lw,label=lab)
    ax[0].set_xscale("log"); ax[0].set_xlabel("tetrahedron volume ÷ mean tetrahedron volume (log scale)",fontsize=FS)
    ax[1].set_xlabel("connections per node (Delaunay degree)",fontsize=FS); ax[2].set_xlabel("link length ÷ mean node spacing",fontsize=FS)
    ttl=["(a) tetrahedron volumes","(b) number of connections","(c) link lengths"]
    for a_,t in zip(ax,ttl): a_.set_ylabel("fraction",fontsize=FS); a_.set_title(t,fontsize=FT,loc="left"); a_.tick_params(labelsize=10.5); a_.grid(alpha=0.25)
    fig.suptitle(title,fontsize=FT+1,fontweight="bold")
    h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc="upper center",bbox_to_anchor=(0.5,0.095),ncol=2,frameon=False,fontsize=FL)
    fig.savefig(D+fn,dpi=140,bbox_inches="tight"); plt.close(fig)
hist3([("fluorite | perfect","perfect FCC + tetra-hole voids","k",1.8,False),("web ratio | random","fully random (Poisson)","#e67e22",1.8,False),
       ("web BEFORE","foam s2 web BEFORE edits","0.45",0,True),("web AFTER","foam s2 web AFTER edits","#1f3a93",1.8,False)],
      "fcc_vs_web_hist.png","Perfect lattice vs random vs foam s2 web (3D Delaunay)")
hist3([("fluorite | perfect","perfect FCC","k",1.8,False),("fluorite | jitter 5%","jitter 5% of spacing","#27ae60",1.6,False),
       ("fluorite | jitter 15%","jitter 15%","#2980b9",1.6,False),("fluorite | jitter 30%","jitter 30%","#8e44ad",1.6,False),("fluorite | random","fully random","#e67e22",1.8,False)],
      "ladder_hist_jitter.png","Imperfect crystal 1: thermal jitter (FCC + tetra-hole voids)")
hist3([("fluorite | perfect","perfect FCC","k",1.8,False),("fluorite | defects 3%","point defects 3% each","#27ae60",1.6,False),
       ("fluorite | defects 10%","point defects 10% each","#2980b9",1.6,False),("fluorite | polycrystal","polycrystal (10 grains)","#8e44ad",1.6,False),("fluorite | random","fully random","#e67e22",1.8,False)],
      "ladder_hist_defects.png","Imperfect crystal 2: point defects and grains (FCC + tetra-hole voids)")
# ---------- entropy vs energy
ORDER=["perfect","jitter 5%","jitter 15%","jitter 30%","defects 3%","defects 10%","polycrystal","random"]
def ladder_plot(comp,fn,title,web):
    fig,ax=plt.subplots(figsize=(9,7.6)); fig.subplots_adjust(bottom=0.25)
    pts={n:(R[f"{comp} | {n}"]["S_hist_logvol"],R[f"{comp} | {n}"]["u_dimless"]) for n in ORDER}
    col={"perfect":"k","jitter 5%":"#27ae60","jitter 15%":"#27ae60","jitter 30%":"#27ae60","defects 3%":"#2980b9","defects 10%":"#2980b9","polycrystal":"#8e44ad","random":"#e67e22"}
    mk={"perfect":"s","jitter 5%":"o","jitter 15%":"o","jitter 30%":"o","defects 3%":"D","defects 10%":"D","polycrystal":"^","random":"X"}
    j=[pts[n] for n in ("perfect","jitter 5%","jitter 15%","jitter 30%","random")]; ax.plot(*np.array(j).T,color="0.75",lw=1,zorder=1)
    off={"jitter 15%":(8,-16),"jitter 30%":(-20,-20),"defects 3%":(-10,10),"defects 10%":(8,6),"polycrystal":(8,8),"perfect":(-10,10)}
    for n,(x,y) in pts.items():
        ax.scatter(x,y,s=110,c=col[n],marker=mk[n],zorder=3); ax.annotate(n,(x,y),xytext=off.get(n,(8,6)),textcoords="offset points",fontsize=FL)
    if web:
        for k,c_,o in (("web BEFORE","0.45",(-150,4)),("web AFTER","#1f3a93",(-145,-6))):
            x,y=R[k]["S_hist_logvol"],R[k]["u_dimless"]; ax.scatter(x,y,s=160,c=c_,marker="*",zorder=4); ax.annotate(k.replace("web","foam s2 web"),(x,y),xytext=o,textcoords="offset points",fontsize=FL,color=c_)
    ax.set_xlabel("disorder entropy: spread of tetrahedron sizes (nats)\n← ordered          disordered →",fontsize=FS)
    ax.set_ylabel("link energy per node, U/N ÷ (G·M²_cluster / spacing)\n← more bound          less bound →",fontsize=FS)
    ax.margins(x=0.08,y=0.1); ax.grid(alpha=0.3); ax.tick_params(labelsize=10.5); ax.set_title(title,fontsize=FT,fontweight="bold")
    fig.text(0.5,0.03,"Entropy = Shannon entropy of the histogram of log(tetrahedron volume), 0.1-dex bins.\nEnergy = −G Σ MᵢMⱼ / rᵢⱼ over Delaunay links, voids as negative mass, per node, in units set by each structure's own spacing.\nGrey line: perfect → jitter 5, 15, 30% → random.",ha="center",fontsize=10.5,color="0.25")
    fig.savefig(D+fn,dpi=140,bbox_inches="tight"); plt.close(fig)
ladder_plot("web ratio","entropy_energy_ladder.png","From perfect crystal to random, with the foam s2 web\n(lattices with the web's mix: 1 void per 6 clusters)",True)
ladder_plot("fluorite","entropy_energy_ladder_fluorite.png","Same ladder with voids in every tetrahedral hole\n(2 voids per cluster; strongly bound by void–void links)",False)
# ---------- polycrystal slice
Z=np.load(D+"polycrystal_fluorite.npz"); P,K,Gr,L=Z["P"],Z["K"],Z["G"],float(Z["L"]); m=np.abs(P[:,2]-L/2)<0.04
_pal=["#f39c12","#27ae60","#8e44ad","#16a085","#d35400","#f1c40f","#b8860b","#e84393","#c71585","#a3cb38"]; cm=lambda i:_pal[i%10]; fig,ax=plt.subplots(figsize=(8.4,9.2)); fig.subplots_adjust(bottom=0.2)
for g in np.unique(Gr[m]):
    s=m&(Gr==g); c_="0.7" if g<0 else cm(g%10)
    ax.scatter(P[s&(K==0),0],P[s&(K==0),1],s=9,c=[c_],lw=0); ax.scatter(P[s&(K==1),0],P[s&(K==1),1],s=9,facecolors="none",edgecolors=[c_],lw=0.6)
ax.set_aspect("equal"); ax.set_xlim(0,L); ax.set_ylim(0,L); ax.set_xlabel("x (Gpc)",fontsize=FS); ax.set_ylabel("y (Gpc)",fontsize=FS); ax.tick_params(labelsize=10.5)
ax.set_title(f"Polycrystal side slice (80 Mpc thick) of a {L:.2f} Gpc box",fontsize=FT,fontweight="bold")
fig.legend(handles=[Patch(color=cm(0),label="each colour = one crystal grain (own orientation, spacing ±5%)"),Patch(color="0.7",label="grey = disordered grain boundary"),
    Line2D([],[],marker="o",ls="",ms=5,color="0.2",label="filled dot = cluster"),Line2D([],[],marker="o",ls="",ms=5,mfc="none",mec="0.2",label="open circle = void")],
    loc="lower center",bbox_to_anchor=(0.5,0.02),ncol=1,frameon=False,fontsize=FL)
fig.savefig(D+"polycrystal_slice.png",dpi=140,bbox_inches="tight"); plt.close(fig)
print("ok")
