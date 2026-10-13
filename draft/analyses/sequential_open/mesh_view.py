"""Mesh view (Chris, 10 Oct): non-crossing network of the side-view slab (|z|<1.5 Gpc), BEFORE and AFTER, for one minimal-edit run.
Network rule = void_mass_network/planar.py, applied in the in-plane projection (so the drawn mesh has no crossings at all):
  candidates: each slab node's 10 nearest slab neighbours in 3D (<= 6 Gpc), typed cluster-cluster / cluster-void / void-void;
  greedy shortest-first; a link is dropped if it passes through another node's body (void slab cross-section, 0.95 R) or if it comes
  within 2 r_t of an already-kept link that shares no end node (r_t = 0.05 x median void radius, as planar.py).
usage: python3 mesh_view.py <minedit_ent json> <out prefix>"""
import sys,json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection,PatchCollection
from matplotlib.patches import Circle,Rectangle
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
sys.path.insert(0,"/workspace/clump_scan/void_mass_network"); import planar as PL
import edit_before_after as EB
from common import *
def network(X,r,v):
    sl=np.where((np.abs(X[:,2])<1.5)&(r>0)&(np.linalg.norm(X,axis=1)<=30))[0]; P=X[sl]; vv=v[sl]
    dd,jj=cKDTree(P).query(P,k=11); a=np.repeat(np.arange(len(P)),10); b=jj[:,1:].ravel(); L3=dd[:,1:].ravel()
    m=(L3<=6.0); a,b=a[m],b[m]; E=np.unique(np.sort(np.c_[a,b],axis=1),axis=0)
    Q=np.c_[P[:,:2],np.zeros(len(P))]; P0=Q[E[:,0]].copy(); P1=Q[E[:,1]].copy(); L=np.linalg.norm(P1-P0,axis=1); mid=0.5*(P0+P1)
    BR=np.where(vv,0.95*np.sqrt(np.maximum(r[sl]**2-P[:,2]**2,0)),0.0)
    bl=cKDTree(Q).query_ball_point(mid,L/2+BR.max()); bptr,bidx=PL._csr([np.array(x,int) for x in bl])
    okb=PL.body_check(P0,P1,E[:,0].astype(np.int64),E[:,1].astype(np.int64),bptr,bidx,Q,BR)
    thr=2*0.05*np.median(r[v&(r>0)])
    nl=cKDTree(mid).query_ball_point(mid,L+thr); nptr,nidx=PL._csr(nl)
    kept,why=PL.greedy(P0,P1,E[:,0].astype(np.int64),E[:,1].astype(np.int64),np.argsort(L,kind="stable"),nptr,nidx,okb,thr)
    E=E[kept]; T=np.where(vv[E[:,0]]&vv[E[:,1]],"vv",np.where(~vv[E[:,0]]&~vv[E[:,1]],"cc","cv"))
    return sl,E,T
def draw(ax,X,r,v,sl,E,hl,title,lim=None):
    cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=0); ax.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=ax.transData))
    G=hl[sl]; mv=G[E[:,0]]|G[E[:,1]]; P=X[sl,:2]
    ax.add_collection(LineCollection(np.stack([P[E[~mv,0]],P[E[~mv,1]]],1),colors="0.55",lw=0.45,zorder=2))
    ax.add_collection(LineCollection(np.stack([P[E[mv,0]],P[E[mv,1]]],1),colors=EB.HL,lw=1.0,zorder=3))
    vs=v[sl]; s=7 if lim is None else 22
    for m_,c in ((~vs&~G,"#c0392b"),(vs&~G,"#2e86c1")): ax.scatter(P[m_,0],P[m_,1],s=s*(1.6 if c=="#2e86c1" else 1),c=c,lw=0,zorder=4)
    ax.scatter(P[G,0],P[G,1],s=s*1.5,c=EB.HL,edgecolor=EB.HLE,lw=0.4,zorder=5)
    th=np.linspace(0,2*np.pi,400); ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=2.0,zorder=6)); ax.plot(0,0,"k+",ms=11,mew=2,zorder=7)
    if lim is None:
        frame30(ax,lw=1.2)
        for t in ax.texts: t.set_fontsize(EB.FK)
        ax.set_xticks([-30,0,30]); ax.set_yticks([-30,0,30])
    else:
        ax.plot(30*np.cos(th),30*np.sin(th),color="0.4",ls=":",lw=1); ax.set_xlim(lim[0]); ax.set_ylim(lim[1]); ax.set_aspect("equal")
        ax.set_xlabel("x (Gpc)",fontsize=EB.FS); ax.set_ylabel("y (Gpc)",fontsize=EB.FS)
    ax.tick_params(labelsize=EB.FK); ax.set_title(title,fontsize=EB.FT)
def legend(fig,zoom=False):
    hd=[Line2D([],[],color="0.55",lw=1,label="link (non-crossing network)"),Line2D([],[],color=EB.HL,lw=1.6,label="link touching a moved body"),
        Line2D([],[],marker="o",ls="",ms=5,color="#c0392b",label="cluster"),Line2D([],[],marker="o",ls="",ms=6,color="#2e86c1",label="void centre"),
        Line2D([],[],marker="o",ls="",ms=6,mfc=EB.HL,mec=EB.HLE,label="moved body (at old spot BEFORE, new spot AFTER)"),
        Line2D([],[],color="#1f4e79",lw=2,label="CMB ring / our sphere")]
    if not zoom: hd+= [Rectangle((0,0),1,1,fc=EB.BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)"),Rectangle((0,0),1,1,fc="#fadbd8",label="CoE core cap (inside 30 Gpc)"),Line2D([],[],color="0.4",ls=":",lw=1.2,label="30 Gpc edge of model"),Line2D([],[],color="#a93226",marker=">",label="direction to the CoE")]
    fig.legend(handles=hd,loc="upper center",ncol=3,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.02),handlelength=1.6,columnspacing=1.3,
               title="nodes and links in the slab |z| < 1.5 Gpc (side-view slab); network built separately for each view",title_fontsize=EB.FK)
if __name__=="__main__":
    js,pre=sys.argv[1],sys.argv[2]; d=json.load(open(js)); base=d["tag"].replace("_entropy","")
    Xb,rb,vb,_,_=EB.layout(base); Xa,ra,va,_,_=EB.layout(base,d)
    idx=np.array(d["idx"],int); HLm=np.zeros(len(Xb),bool); HLm[idx[np.array(d["moved"],bool)]]=True
    nets=[network(Xb,rb,vb),network(Xa,ra,va)]; info={}
    for nm,(sl,E,T) in zip(("before","after"),nets):
        G=HLm[sl]; info[nm]=dict(nodes=len(sl),links=len(E),cc=int((T=="cc").sum()),cv=int((T=="cv").sum()),vv=int((T=="vv").sum()),moved_nodes=int(G.sum()),links_moved=int((G[E[:,0]]|G[E[:,1]]).sum()))
    print(json.dumps(info)); json.dump(info,open(pre+"_links.json","w"),indent=1)
    lab=lambda k: f"{k.upper()}: {info[k]['links']:,} links, {info[k]['nodes']:,} nodes\n{info[k]['links_moved']:,} links touch {info[k]['moved_nodes']} moved bodies"
    fig=plt.figure(figsize=(13,10.4)); GS=fig.add_gridspec(2,2,height_ratios=[1.0,0.42],hspace=0.36,wspace=0.1,left=0.06,right=0.98,top=0.9,bottom=0.08); ax=[fig.add_subplot(GS[0,0]),fig.add_subplot(GS[0,1])]
    EB.residual_below(fig,GS[1,:],d,EB.layout(base)[4])
    fig.suptitle("foam growth, seed 2 — entropy-rule edits as a mesh (non-crossing network)",fontsize=EB.FT+1,fontweight="bold",color="#8e44ad")
    for a_,(X,r,v),(sl,E,T),k in zip(ax,((Xb,rb,vb),(Xa,ra,va)),nets,("before","after")): draw(a_,X,r,v,sl,E,HLm,lab(k))
    legend(fig); fig.savefig(pre+"_before_after.png",dpi=140,bbox_inches="tight"); plt.close(fig)
    # zoom: 10x10 Gpc window with the most moved slab nodes, centred near the ring
    Pm=Xa[HLm&(np.abs(Xa[:,2])<1.5)]; best=None
    for t_ in np.linspace(0,2*np.pi,72,endpoint=False):
        c=RMIN*np.array([np.cos(t_),np.sin(t_)]); nI=np.sum(np.all(np.abs(Pm[:,:2]-c)<5,axis=1))
        if best is None or nI>best[0]: best=(nI,c)
    c=best[1]; lim=((c[0]-5,c[0]+5),(c[1]-5,c[1]+5))
    fig,ax=plt.subplots(1,2,figsize=(13,7.0)); fig.subplots_adjust(wspace=0.22,left=0.07,right=0.98,top=0.86,bottom=0.12)
    fig.suptitle(f"Zoom: busiest moved region on the ring (10 × 10 Gpc around x={c[0]:.1f}, y={c[1]:.1f} Gpc)",fontsize=EB.FT+1,fontweight="bold",color="#8e44ad")
    for a_,(X,r,v),(sl,E,T),k in zip(ax,((Xb,rb,vb),(Xa,ra,va)),nets,("before","after")): draw(a_,X,r,v,sl,E,HLm,k.upper(),lim=lim)
    legend(fig,zoom=True); fig.savefig(pre+"_zoom.png",dpi=140,bbox_inches="tight"); print("zoom",best)
