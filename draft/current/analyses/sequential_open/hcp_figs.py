"""Figures for the HCP lattice analogue (Chris, 11 Oct): hcp_cell.png, hcp_vs_fcc.png, hcp_slab.png, hcp_entropy_energy_ladder.png,
hcp_vs_fcc_hist.png, hcp_shells.png. Run after hcp_lattice.py / hcp_seeds.py / hcp_ewald.py. Line style: cluster-cluster links SOLID dark grey,
cluster-void links DASHED blue; clusters red, voids blue. Moderate fonts, legends below."""
import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import cKDTree
D="hcp_lattice/"; FS=12; FT=13; FL=11
RC,REC,VC,VEC,LK,VL="#c0392b","#7b241c","#2e86c1","#1b4f72","0.25","#2e86c1"
d=171/np.sqrt(2); c=np.sqrt(8/3)*d; h=np.sqrt(3)/2*d
def tri_layer(off,R):
    """triangular-lattice positions (spacing d) shifted by off (2-vector), within radius R of the axis"""
    n=int(R/d)+3; P=np.array([[i*d+j*d/2,j*h] for i in range(-n,n+1) for j in range(-n,n+1)])+off
    return P[np.hypot(P[:,0],P[:,1])<=R+1e-6]
OFF={"A":np.array([0,0]),"B":np.array([d/2,d/(2*np.sqrt(3))]),"C":np.array([0,d/np.sqrt(3)])}
def stack(seq,R):
    P=[];L=[]
    for k,s in enumerate(seq):
        q=tri_layer(OFF[s],R); P.append(np.c_[q,np.full(len(q),k*c/2)]); L+= [k]*len(q)
    return np.vstack(P),np.array(L)
def nn_pairs(P,r=d):
    t=cKDTree(P); return np.array(sorted(t.query_pairs(r*1.01)))
def holes(P):
    """tetra/octa holes of a close-packed stack: centres of empty spheres touching 4 (r=0.612d) / 6 (r=0.707d) clusters"""
    t=cKDTree(P); out={"tetra":[],"octa":[]}
    for k in range(int(P[:,2].max()/(c/2))):
        z0=k*c/2
        for zf,nm,rr,n in ((0.25,"octa",np.sqrt(.5)*d,6),(0.375,"tetra",np.sqrt(3/8)*d,4),(0.625,"tetra",np.sqrt(3/8)*d,4)):
            pass
    # brute force candidate positions: all A/B/C columns at the hole heights
    cand=[]
    for k in range(int(round(P[:,2].max()/(c/2)))):
        for s in "ABC":
            q=tri_layer(OFF[s],P[:,0].max()+d)
            for f in (0.25,0.5,0.75): cand.append(np.c_[q,np.full(len(q),(k+f)*c/2)])
    cand=np.vstack(cand)
    for nm,rr,n in (("tetra",np.sqrt(3/8)*d,4),("octa",np.sqrt(.5)*d,6)):
        dd,ii=t.query(cand,n+1)
        ok=(np.abs(dd[:,:n]-rr)<1e-3*d).all(1)&(dd[:,n]>rr+1e-2*d)
        out[nm]=[(cand[j],ii[j,:n]) for j in np.where(ok)[0]]
    return out
# ======================= 1) hcp_cell.png
R1=d*1.001; P,Ly=stack("ABA",R1)
# keep the hexagonal prism: A layers = centre + 6 corners, B layer = the 3 inner sites
inP=np.hypot(P[:,0],P[:,1])<=d*1.001; P,Ly=P[inP],Ly[inP]
extra=np.array([[0,2*d/np.sqrt(3),c/2]])   # B cluster of the neighbouring prism, needed to close the octahedron
Pall=np.vstack([P,extra])
fig=plt.figure(figsize=(19,8.6)); gs=fig.add_gridspec(1,3,width_ratios=[1,1,0.95],wspace=0.14,left=0.01,right=0.99,top=0.88,bottom=0.22)
E=nn_pairs(Pall); th=np.deg2rad(np.arange(0,361,60))
def prism(ax,which):
    for i,j in E: ax.plot(*np.c_[Pall[i],Pall[j]],color=LK,lw=1.5,zorder=2,alpha=1 if (which=="octa" or not ((Pall[i]==extra[0]).all() or (Pall[j]==extra[0]).all())) else 0)
    for t_ in th[:-1]: ax.plot([d*np.cos(t_)]*2,[d*np.sin(t_)]*2,[0,c],color="0.7",lw=0.9,ls=":")
    if which=="tetra":
        Q=Pall[tn]; ax.add_collection3d(Poly3DCollection([Q[[0,1,2]],Q[[0,1,3]],Q[[0,2,3]],Q[[1,2,3]]],fc="#aed6f1",alpha=0.35,ec="none"))
        for q in Q: ax.plot(*np.c_[tv,q],color=VL,lw=2.0,ls=(0,(4,3)))
        ax.scatter(*P.T,s=240,c=RC,edgecolor=REC,depthshade=False); ax.scatter(*Q.T,s=240,c=RC,edgecolor="k",linewidths=2,depthshade=False)
        ax.scatter(*tv,s=330,c=VC,edgecolor=VEC,depthshade=False)
        ax.set_title("(a) void in a tetrahedral hole\n4 clusters: 1 below (layer A) + 3 above (layer B)",fontsize=FT)
    else:
        ax.add_collection3d(Poly3DCollection(f_,fc="#aed6f1",alpha=0.3,ec="none"))
        for q in Qo: ax.plot(*np.c_[ov,q],color=VL,lw=2.0,ls=(0,(4,3)))
        ax.scatter(*P.T,s=240,c=RC,edgecolor=REC,depthshade=False); ax.scatter(*extra.T,s=240,c="#f1948a",edgecolor=REC,depthshade=False)
        ax.scatter(*Qo.T,s=240,facecolors="none",edgecolor="k",linewidths=2,depthshade=False)
        ax.scatter(*ov,s=330,c=VC,edgecolor=VEC,depthshade=False)
        ax.set_title("(b) void in an octahedral hole\n6 clusters: 3 in layer A + 3 in layer B",fontsize=FT)
    for z,s_ in ((0,"A"),(c/2,"B"),(c,"A")): ax.text(1.3*d,-0.2*d,z,f"layer {s_}",fontsize=FL,fontweight="bold",color="0.2")
    ax.set_xlabel("x (Mpc)",fontsize=FL); ax.set_ylabel("y (Mpc)",fontsize=FL); ax.set_zlabel("z (Mpc)",fontsize=FL,labelpad=4)
    ax.set_xticks([-121,0,121]); ax.set_yticks([-121,0,121]); ax.set_zticks([0,99,198]); ax.tick_params(labelsize=9.5)
    ax.view_init(elev=13,azim=-62); ax.set_proj_type("ortho"); ax.set_box_aspect((np.ptp(Pall[:,0]),np.ptp(Pall[:,1]),c),zoom=1.12)
tv=np.array([0,0,3*c/8]); HL=holes(Pall)
ti=[k for k,(x,_) in enumerate(HL["tetra"]) if np.linalg.norm(x-tv)<1][0]; tv,tn=HL["tetra"][ti]
ov=np.array([0,d/np.sqrt(3),c/4]); oi=[k for k,(x,_) in enumerate(HL["octa"]) if np.linalg.norm(x-ov)<1][0]; ov,on=HL["octa"][oi]
Qo=Pall[on]; f_=[]
for a_ in range(6):
    for b_ in range(a_+1,6):
        for e_ in range(b_+1,6):
            if all(abs(np.linalg.norm(Qo[x]-Qo[y])-d)<1 for x,y in ((a_,b_),(a_,e_),(b_,e_))): f_.append(Qo[[a_,b_,e_]])
prism(fig.add_subplot(gs[0],projection="3d"),"tetra"); prism(fig.add_subplot(gs[1],projection="3d"),"octa")
fig.suptitle("One HCP unit: hexagonal prism, layers A-B-A; neighbour spacing d = 121 Mpc, layer repeat c = 198 Mpc",fontsize=FT+1,fontweight="bold")
# top view
a2=fig.add_subplot(gs[2]); Rt=2.2*d
for s,kw in (("A",dict(s=330,c=RC,edgecolor=REC,zorder=3)),("B",dict(s=200,facecolors="#f5b7b1",edgecolors=REC,linewidths=1.5,zorder=4))):
    q=tri_layer(OFF[s],Rt); a2.scatter(q[:,0],q[:,1],**kw)
q=tri_layer(OFF["C"],Rt); a2.scatter(q[:,0],q[:,1],s=60,marker="x",c="0.45",zorder=4,lw=1.5)
qa=tri_layer(OFF["A"],Rt); 
for i,j in nn_pairs(np.c_[qa,np.zeros(len(qa))]): a2.plot(*np.c_[qa[i],qa[j]],color="0.75",lw=0.8,zorder=1)
hx=np.c_[d*np.cos(th),d*np.sin(th)]; a2.plot(hx[:,0],hx[:,1],color="0.2",lw=1.2,ls="--",zorder=5)
for s,dy in (("A",-0.3),("B",-0.3),("C",-0.3)):
    p=OFF[s]+([d,0] if s!="C" else [d,0]); a2.annotate(s,p,xytext=(0,-22),textcoords="offset points",ha="center",fontsize=FS,fontweight="bold")
a2.yaxis.tick_right(); a2.yaxis.set_label_position("right"); a2.set_aspect("equal"); a2.set_xlim(-Rt,Rt); a2.set_ylim(-Rt,Rt); a2.set_xlabel("x (Mpc)",fontsize=FL); a2.set_ylabel("y (Mpc)",fontsize=FL); a2.tick_params(labelsize=10)
a2.set_title("(c) top view of the stacking: A B A B …\nC sites stay empty in HCP (FCC fills them: A B C …)",fontsize=FT)
fig.legend(handles=[Line2D([],[],marker="o",ls="",ms=12,mfc=RC,mec=REC,label="cluster (black ring = the void's own clusters)"),
    Line2D([],[],marker="o",ls="",ms=12,mfc="#f1948a",mec=REC,label="B cluster of the next prism (closes the octahedron)"),
    Line2D([],[],marker="o",ls="",ms=12,mfc=VC,mec=VEC,label="void in a tetrahedral hole (r = 0.61 d = 74 Mpc to its clusters)"),
    Line2D([],[],marker="o",ls="",ms=12,mfc=VC,mec=VEC,label="void in an octahedral hole (r = 0.71 d = 85 Mpc to its clusters)"),
    Line2D([],[],color=LK,lw=1.8,label="cluster–cluster link (solid, length d)"),
    Line2D([],[],color=VL,lw=1.8,ls=(0,(4,3)),label="cluster–void link (dashed)"),
    Line2D([],[],marker="o",ls="",ms=11,mfc=RC,mec=REC,label="top view: layer A"),Line2D([],[],marker="o",ls="",ms=9,mfc="#f5b7b1",mec=REC,label="top view: layer B"),
    Line2D([],[],marker="x",ls="",ms=8,color="0.45",mew=1.5,label="top view: empty C sites (octahedral holes stack above them)"),
    Line2D([],[],color="0.2",ls="--",lw=1.2,label="top view: outline of the prism")],
    loc="lower center",ncol=3,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.0),columnspacing=2,
    title="per cluster: 2 tetrahedral holes and 1 octahedral hole (same as FCC); one of each is filled here",title_fontsize=FL)
fig.savefig(D+"hcp_cell.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# ======================= 2) hcp_vs_fcc.png
def stack_panel(ax3,ax2,seq,title):
    Rs=1.5*d; P,Ly=stack(seq,Rs); E=nn_pairs(P); H=holes(P)["octa"]
    for i,j in E: ax3.plot(*np.c_[P[i],P[j]],color=LK,lw=1.0)
    for x,nb in H:
        for q in P[nb]: ax3.plot(*np.c_[x,q],color=VL,lw=0.9,ls=(0,(3,2.5)))
    ax3.scatter(*P.T,s=110,c=RC,edgecolor=REC,depthshade=False)
    if H: Xh=np.array([x for x,_ in H]); ax3.scatter(*Xh.T,s=120,c=VC,edgecolor=VEC,depthshade=False)
    for k,s in enumerate(seq): ax3.text(Rs*1.15,-Rs*0.6,k*c/2,s,fontsize=FS,fontweight="bold")
    ax3.view_init(elev=12,azim=-65); ax3.set_proj_type("ortho"); ax3.set_box_aspect((2*Rs,2*Rs,c/2*(len(seq)-1)))
    ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks(np.arange(len(seq))*c/2); ax3.set_zticklabels([f"{z:.0f}" for z in np.arange(len(seq))*c/2],fontsize=9)
    ax3.set_zlabel("height (Mpc)",fontsize=FL); ax3.set_title(title,fontsize=FT,fontweight="bold")
    # side view: thin slab through the axis, looking along x (A, B, C columns separate in y)
    Rw=2.6*d; P2,_=stack(seq,Rw); P2=P2[np.abs(P2[:,0])<0.3*d]; E2=nn_pairs(P2)
    Hs=holes(stack(seq,Rw+d)[0]); Pf=stack(seq,Rw+d)[0]
    for i,j in E2: ax2.plot(*np.c_[P2[i,1:],P2[j,1:]],color=LK,lw=1.2,zorder=2)
    Ho=[(x,nb) for x,nb in Hs["octa"] if abs(x[0])<0.3*d and abs(x[1])<Rw-0.8*d]
    for x,nb in Ho:
        for q in Pf[nb]:
            if abs(q[0])<0.3*d: ax2.plot([x[1],q[1]],[x[2],q[2]],color=VL,lw=1.0,ls=(0,(3,2.5)),zorder=1)
    ax2.scatter(P2[:,1],P2[:,2],s=90,c=RC,edgecolor=REC,zorder=3)
    if Ho: Xh=np.array([x for x,_ in Ho]); ax2.scatter(Xh[:,1],Xh[:,2],s=95,c=VC,edgecolor=VEC,zorder=4)
    for k,s in enumerate(seq): ax2.text(Rw+0.15*d,k*c/2,s,fontsize=FS,fontweight="bold",va="center")
    ax2.set_aspect("equal"); ax2.set_xlim(-Rw-0.1*d,Rw+0.5*d); ax2.set_ylim(-0.4*c,(len(seq)-0.6)*c/2)
    ax2.set_xlabel("y (Mpc)",fontsize=FL); ax2.set_ylabel("height (Mpc)",fontsize=FL); ax2.tick_params(labelsize=10)
fig=plt.figure(figsize=(14,12.5)); g=fig.add_gridspec(2,2,height_ratios=[1.25,1],hspace=0.18,wspace=0.12,left=0.06,right=0.97,top=0.93,bottom=0.17)
stack_panel(fig.add_subplot(g[0,0],projection="3d"),fig.add_subplot(g[1,0]),"ABABA","HCP: A B A B (every 2nd layer repeats)")
stack_panel(fig.add_subplot(g[0,1],projection="3d"),fig.add_subplot(g[1,1]),"ABCAB","FCC: A B C A B C (every 3rd layer repeats)")
fig.text(0.5,0.565,"Side views: a thin slice through the middle, looking along x",ha="center",fontsize=FT-1,style="italic")
fig.suptitle("HCP vs FCC at the same spacing (d = 121 Mpc), voids in all octahedral holes",fontsize=FT+1,fontweight="bold")
fig.legend(handles=[Line2D([],[],marker="o",ls="",ms=10,mfc=RC,mec=REC,label="cluster"),Line2D([],[],marker="o",ls="",ms=10,mfc=VC,mec=VEC,label="void (octahedral hole, 6 clusters)"),
    Line2D([],[],color=LK,lw=1.5,label="cluster–cluster link (solid)"),Line2D([],[],color=VL,lw=1.5,ls=(0,(3,2.5)),label="cluster–void link (dashed)")],
    loc="lower center",ncol=4,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.085))
fig.text(0.5,0.02,"The first two layers are identical in both. The difference is the third layer: HCP puts it straight above the first (A), FCC on the empty C sites.\n"
         "So in HCP the octahedral holes (voids) line up in straight vertical columns, sharing faces, 0.82 d apart;\nin FCC they step sideways one site per layer (a staircase) and are never closer than 1.0 d.",ha="center",fontsize=FL,color="0.25")
fig.savefig(D+"hcp_vs_fcc.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("figs 1-2 ok")
# ======================= 3) hcp_slab.png : 1 Gpc side view
exec("HB="+open("hcp_lattice.py").read().split("HB=")[1].split("def crystal_hcp")[0])   # hex_cart, hcp_sites (Gpc units)
dG=d/1e3; rng=np.random.default_rng(5); W=1.0
cs,ts,os_=hcp_sites(np.array([-.05,-1.5*dG,-.05]),np.array([W+.05,1.5*dG,W+.05]),dG)
fig,axs=plt.subplots(1,2,figsize=(15,8.6)); fig.subplots_adjust(wspace=0.12,left=0.06,right=0.98,top=0.86,bottom=0.2)
for ax,(H,fill,nm,k) in zip(axs,((ts,1/12,"voids in 1/12 of the tetrahedral holes",4),(os_,1/6,"voids in 1/6 of the octahedral holes",6))):
    sl=lambda X,w=0.5: np.abs(X[:,1])<w*dG
    C=cs[sl(cs)]; Hv=H[sl(H)&(rng.random(len(H))<fill)]; Hv=Hv[np.all((Hv[:,[0,2]]>=0)&(Hv[:,[0,2]]<=W),axis=1)]
    t=cKDTree(C)
    for i,j in t.query_pairs(dG*1.01): ax.plot([C[i,0],C[j,0]],[C[i,2],C[j,2]],color=LK,lw=0.8,zorder=2)
    T=cKDTree(cs); dd,ii=T.query(Hv,k); out=set()
    for x,nb in zip(Hv,ii):
        for q in nb:
            ax.plot([x[0],cs[q,0]],[x[2],cs[q,2]],color=VL,lw=1.0,ls=(0,(3,2)),zorder=3)
            if not sl(cs[q:q+1])[0]: out.add(q)
    m=np.all((C[:,[0,2]]>=-0.01)&(C[:,[0,2]]<=W+0.01),axis=1)
    ax.scatter(C[m,0],C[m,2],s=22,c=RC,lw=0,zorder=4)
    if out: O=cs[list(out)]; ax.scatter(O[:,0],O[:,2],s=22,facecolors="none",edgecolors=RC,lw=0.9,zorder=4)
    ax.scatter(Hv[:,0],Hv[:,2],s=70,c=VC,edgecolor=VEC,lw=0.6,zorder=5)
    ax.set_aspect("equal"); ax.set_xlim(0,W); ax.set_ylim(0,W); ax.set_xlabel("x (Gpc)",fontsize=FS); ax.set_ylabel("height along the stacking axis (Gpc)",fontsize=FS)
    ax.tick_params(labelsize=10.5); ax.set_title(f"{nm}\n{m.sum()} clusters, {len(Hv)} voids in view (1 void per 6 clusters overall)",fontsize=FT)
fig.suptitle("HCP block, 1 Gpc wide, side view: slab 121 Mpc thick (one A row + one B row)",fontsize=FT+1,fontweight="bold")
fig.legend(handles=[Line2D([],[],marker="o",ls="",ms=6,mfc=RC,mec=RC,label="cluster in the slab"),Line2D([],[],marker="o",ls="",ms=6,mfc="none",mec=RC,label="cluster just outside the slab (a void's own cluster)"),
    Line2D([],[],marker="o",ls="",ms=9,mfc=VC,mec=VEC,label="void (lattice hole)"),Line2D([],[],color=LK,lw=1.2,label="cluster–cluster link (solid)"),Line2D([],[],color=VL,lw=1.2,ls=(0,(3,2)),label="cluster–void link (dashed)")],
    loc="lower center",ncol=3,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.04))
fig.text(0.5,0.0,"Lattice spacing is 121 Mpc; the web's mean spacing is 1.9 Gpc, so this whole block would fit inside one web cell. Only the pattern compares, not the scale.",ha="center",fontsize=FL-0.5,color="0.3")
fig.savefig(D+"hcp_slab.png",dpi=130,bbox_inches="tight"); plt.close(fig)
# ======================= 4) ladder
RF=json.load(open("fcc_lattice/ladder.json")); RH=json.load(open(D+"hcp_ladder.json")); SD=json.load(open(D+"hcp_fill_seeds.json"))
P_=lambda o:(o["S_hist_logvol"],o["u_dimless"])
fig,ax=plt.subplots(figsize=(11,8.4)); fig.subplots_adjust(bottom=0.3)
names=["perfect","jitter 5%","jitter 15%","jitter 30%","defects 3%","defects 10%","polycrystal"]
mk={"perfect":"s","jitter 5%":"o","jitter 15%":"o","jitter 30%":"o","defects 3%":"D","defects 10%":"D","polycrystal":"^"}
F={n:P_(RF[f"web ratio | {n}"]) for n in names}; Hh={n:P_(RH[f"HCP web ratio | {n}"]) for n in names}; rnd=P_(RF["web ratio | random"])
for Q,col,lw in ((F,"0.7",1),(Hh,"#c0392b",1.2)):
    ax.plot(*np.array([Q[n] for n in ("perfect","jitter 5%","jitter 15%","jitter 30%")]+[rnd]).T,color=col,lw=lw,zorder=1,ls="-" if col!="0.7" else "--")
for n in names:
    ax.scatter(*F[n],s=90,marker=mk[n],facecolors="none",edgecolors="0.45",lw=1.4,zorder=2)
    ax.scatter(*Hh[n],s=110,marker=mk[n],c="#c0392b",edgecolor="#7b241c",zorder=3)
off={"perfect":(-14,12),"jitter 5%":(6,8),"jitter 15%":(-62,-17),"jitter 30%":(10,-4),"defects 3%":(-30,12),"defects 10%":(-80,8),"polycrystal":(8,8)}
for n in names: ax.annotate(n,Hh[n],xytext=off[n],textcoords="offset points",fontsize=FL-0.5,color="#7b241c")
# fill-seed scatter of the perfect lattices
for k,col in (("HCP tetra 1/12","#c0392b"),("FCC tetra 1/12","0.45")):
    ax.errorbar(SD[k]["S_hist_mean"],SD[k]["u_mean"],yerr=SD[k]["u_sd"],fmt="none",ecolor=col,elinewidth=2,capsize=5,zorder=2)
ho=P_(RH["HCP web ratio octa | perfect"]); ax.scatter(*ho,s=110,marker="P",c="#e74c3c",edgecolor="#7b241c",zorder=3); ax.annotate("HCP, voids in octa holes",ho,xytext=(8,-4),textcoords="offset points",fontsize=FL-0.5,color="#7b241c")
ax.scatter(*rnd,s=120,marker="X",c="#e67e22",zorder=3); ax.annotate("random",rnd,xytext=(8,4),textcoords="offset points",fontsize=FL-0.5)
for k,c_,o,lab in (("web BEFORE","0.45",(-12,10),"foam s2 web BEFORE"),("web AFTER","#1f3a93",(-12,-18),"foam s2 web AFTER")):
    x,y=P_(RF[k]); ax.scatter(x,y,s=200,c=c_,marker="*",zorder=4); ax.annotate(lab,(x,y),xytext=o,textcoords="offset points",fontsize=FL-0.5,color=c_,ha="right")
for k,c_,o,lab in (("comp web s3 BEFORE","#7d3c98",(10,6),"balanced web s3 BEFORE (789×)"),("comp web s3 AFTER","#117a65",(10,-14),"balanced web s3 AFTER (0.998×)")):
    x,y=P_(RH[k]); ax.scatter(x,y,s=200,c=c_,marker="*",zorder=4); ax.annotate(lab,(x,y),xytext=o,textcoords="offset points",fontsize=FL-0.5,color=c_)
ax.set_xlabel("disorder entropy: spread of tetrahedron sizes (nats)\n← ordered          disordered →",fontsize=FS)
ax.set_ylabel("link energy per node, U/N ÷ (G·M²_cluster / spacing)\n← more bound          less bound →",fontsize=FS)
ax.margins(x=0.1,y=0.1); ax.grid(alpha=0.3); ax.tick_params(labelsize=10.5)
ax.set_title("HCP added to the crystal-to-random ladder\n(lattices with the web's mix: 1 void per 6 clusters)",fontsize=FT,fontweight="bold")
fig.legend(handles=[Line2D([],[],marker="s",ls="",ms=9,mfc="#c0392b",mec="#7b241c",label="HCP (filled, red line: perfect → jitter → random)"),
    Line2D([],[],marker="s",ls="",ms=9,mfc="none",mec="0.45",mew=1.4,label="FCC, same recipe (open, grey dashed line)"),
    Line2D([],[],marker="o",ls="",ms=8,color="0.3",label="circle = jitter 5 / 15 / 30%"),Line2D([],[],marker="D",ls="",ms=8,color="0.3",label="diamond = defects 3 / 10%"),
    Line2D([],[],marker="^",ls="",ms=8,color="0.3",label="triangle = polycrystal (10 grains)"),Line2D([],[],color="0.3",lw=2,label="error bar = ±1 s.d. over 16 random void fills (perfect lattice)"),
    Line2D([],[],marker="*",ls="",ms=13,color="0.3",label="stars = our webs (foam s2, balanced s3)"),Line2D([],[],marker="X",ls="",ms=9,color="#e67e22",label="random (Poisson), shared")],
    loc="lower center",ncol=2,frameon=False,fontsize=FL-0.5,bbox_to_anchor=(0.5,0.075))
fig.text(0.5,0.0,"Entropy = Shannon entropy of the histogram of log(tetrahedron volume), 0.1-dex bins. Energy = −G Σ MᵢMⱼ / rᵢⱼ over Delaunay links,\nvoids as negative mass, per node, in units set by each structure's own spacing.",ha="center",fontsize=10.5,color="0.25")
fig.savefig(D+"hcp_entropy_energy_ladder.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# ======================= 5) histograms HCP vs FCC
RWh=np.load(D+"hcp_ladder_raw.npz"); RWf=np.load("fcc_lattice/ladder_raw.npz")
ser=[(RWf,"FCC clusters only","FCC, clusters only","0.55",True),(RWh,"HCP clusters only","HCP, clusters only","#c0392b",False),
     (RWf,"web ratio | perfect","FCC + voids (web mix)","#5dade2",True),(RWh,"HCP web ratio | perfect","HCP + voids (web mix)","#1b4f72",False)]
ser[0]=(RWh,"FCC clusters only","FCC, clusters only","0.55",True)
fig,ax=plt.subplots(3,1,figsize=(10,12.5)); fig.subplots_adjust(hspace=0.45,bottom=0.12,top=0.93)
vb=np.logspace(-2.5,1,60); lb=np.linspace(0.3,2.2,58); db=np.arange(5.5,30.5,1)
for Z,k,lab,col,fill in ser:
    for a_,q,b in ((ax[0],"vol",vb),(ax[1],"L",lb),(ax[2],"deg",db)):
        x=Z[f"{k}|{q}"]; w=np.ones(len(x))/len(x)
        if fill: a_.hist(x,b,weights=w,color=col,alpha=0.4,label=lab)
        else: a_.hist(x,b,weights=w,histtype="step",color=col,lw=2,label=lab)
ax[0].set_xscale("log"); ax[0].set_xlabel("tetrahedron volume ÷ mean tetrahedron volume (log scale)",fontsize=FS)
ax[1].set_xlabel("link length ÷ mean node spacing",fontsize=FS); ax[2].set_xlabel("connections per node (Delaunay degree)",fontsize=FS)
for a_,t in zip(ax,["(a) tetrahedron volumes","(b) link lengths","(c) connections per node"]): a_.set_ylabel("fraction",fontsize=FS); a_.set_title(t,fontsize=FT,loc="left"); a_.tick_params(labelsize=10.5)
fig.suptitle("HCP vs FCC, perfect lattices, same spacing (3D Delaunay)",fontsize=FT+1,fontweight="bold")
h_,l_=ax[0].get_legend_handles_labels(); fig.legend(h_,l_,loc="lower center",ncol=2,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.02))
fig.savefig(D+"hcp_vs_fcc_hist.png",dpi=130,bbox_inches="tight"); plt.close(fig)
# ======================= 6) neighbour shells
SH=json.load(open(D+"hcp_shells.json"))
fig,ax=plt.subplots(1,2,figsize=(14,6)); fig.subplots_adjust(wspace=0.22,bottom=0.25,top=0.85)
for a_,key,t in ((ax[0],"cluster->cluster","(a) clusters around a cluster"),(ax[1],"tetra void->tetra void","(b) tetrahedral holes around a tetrahedral hole")):
    for nm,col,dx in (("FCC","0.55",-0.012),("HCP","#c0392b",0.012)):
        r=np.array([float(x) for x in SH[nm][key]]); n=np.array(list(SH[nm][key].values()))
        a_.bar(r+dx,n,width=0.022,color=col,label=nm)
        for x,y in zip(r,n): a_.text(x+dx*1.6,y+0.3,f"{y:.0f}",ha="right" if dx<0 else "left",fontsize=9.5,color=col)
    a_.set_xlabel("distance ÷ nearest-neighbour spacing d",fontsize=FS); a_.set_ylabel("number of neighbours",fontsize=FS); a_.set_title(t,fontsize=FT); a_.tick_params(labelsize=10.5)
fig.suptitle("Neighbour shells: identical out to √2 d, different from the third shell on",fontsize=FT+1,fontweight="bold")
fig.legend(*ax[0].get_legend_handles_labels(),loc="lower center",ncol=2,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.05))
fig.savefig(D+"hcp_shells.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("all figs ok")
