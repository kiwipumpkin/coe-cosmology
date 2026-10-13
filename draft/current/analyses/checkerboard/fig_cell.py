import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
LAM=171.4; C=LAM/2; RV=14300.0
fig=plt.figure(figsize=(10,17)); gs=fig.add_gridspec(3,2,height_ratios=[1.1,1.0,1.0],hspace=0.38,wspace=0.15)
# (a) Einasto cell: void at centre, superclusters on corners/edges
ax=fig.add_subplot(gs[0,0],projection="3d")
L=LAM; e=[(0,0,0),(L,0,0),(0,L,0),(0,0,L),(L,L,0),(L,0,L),(0,L,L),(L,L,L)]
for a in e:
    for b in e:
        if sum(abs(np.array(a)-np.array(b))>0)==1 and a<b: ax.plot(*zip(a,b),color="#d62728",lw=3.5,alpha=0.55)
P=np.array(e); ax.scatter(*P.T,s=260,c="#d62728",depthshade=False)
u,v=np.mgrid[0:2*np.pi:30j,0:np.pi:15j]; rv=0.36*L
ax.plot_surface(L/2+rv*np.cos(u)*np.sin(v),L/2+rv*np.sin(u)*np.sin(v),L/2+rv*np.cos(v),color="#3b6fd8",alpha=0.25,lw=0)
rng=np.random.default_rng(1)
for a in e[:4]:
    pass
ax.set_box_aspect((1,1,1)); ax.set_xticks([0,85.7,171.4]); ax.set_yticks([0,85.7,171.4]); ax.set_zticks([0,85.7,171.4])
ax.set_xlabel("Mpc"); ax.set_title("(a) One Einasto cell, edge = 1 period = 171 Mpc\nvoid at centre (blue, radius ≈ 60 Mpc), rich superclusters\nat corners (dots) and chains along edges (red)",fontsize=9.5)
# (b) chessboard cell block 2x2x2
ax=fig.add_subplot(gs[0,1],projection="3d")
def cube(o,s):
    x,y,z=o; v=np.array([[x,y,z],[x+s,y,z],[x+s,y+s,z],[x,y+s,z],[x,y,z+s],[x+s,y,z+s],[x+s,y+s,z+s],[x,y+s,z+s]])
    f=[[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]; return [v[i] for i in f]
for i in range(2):
    for j in range(2):
        for k in range(2):
            col="#d62728" if (i+j+k)%2==0 else "#3b6fd8"
            ax.add_collection3d(Poly3DCollection(cube((i*C,j*C,k*C),C*0.97),facecolor=col,alpha=0.45,edgecolor="k",lw=0.4))
ax.set_xlim(0,LAM); ax.set_ylim(0,LAM); ax.set_zlim(0,LAM); ax.set_box_aspect((1,1,1))
ax.set_xticks([0,85.7,171.4]); ax.set_yticks([0,85.7,171.4]); ax.set_zticks([0,85.7,171.4]); ax.set_xlabel("Mpc")
ax.set_title("(b) Chessboard version: 8 cubes of 85.7 Mpc\n(= ½ period) alternating over (red) / under (blue);\nsame 171 Mpc repeat along each axis",fontsize=9.5)
# (c) length ladder
ax=fig.add_subplot(gs[1,:])
items=[("Abell cluster radius (1.5 h⁻¹)",2.1),("max cluster–cluster link inside a supercluster (24 h⁻¹)",34.3),("void radius in the Einasto cell (~0.36 × step)",62),
       ("chessboard cell edge = ½ period",85.7),("Einasto step / full period (120 h⁻¹)",171.4),("BAO scale, for comparison",147),
       ("Gpc",1000),("our radius to the CMB ring",14300)]
y=np.arange(len(items))[::-1]; ax.barh(y,[v for _,v in items],color=["0.6","0.6","#3b6fd8","#d62728","#d62728","0.75","0.85","k"])
ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels([n for n,_ in items],fontsize=9)
for yy,(n,v) in zip(y,items): ax.text(v*1.1,yy,f"{v:,.0f} Mpc" if v>=10 else f"{v:.1f} Mpc",va="center",fontsize=9)
ax.set_xlim(1,1e5); ax.set_xlabel("length (Mpc, h = 0.7)"); ax.set_title("(c) How a cell scales against clusters, voids and our sphere",fontsize=10)
# (d) mapping to the ring
ax=fig.add_subplot(gs[2,:]); ax.axis("off")
circ=2*np.pi*RV
rows=[["","chessboard cell (85.7 Mpc)","full period (171.4 Mpc)"],
 ["cells / periods from us to the ring",f"{RV/C:.0f}",f"{RV/LAM:.0f}"],
 ["across our sphere (diameter)",f"{2*RV/C:.0f}",f"{2*RV/LAM:.0f}"],
 ["around the ring circumference",f"{circ/C:,.0f}",f"{circ/LAM:,.0f}"],
 ["angular size seen at the ring",f"{np.degrees(C/RV):.2f}°  ({np.degrees(C/RV)*60:.0f}′)",f"{np.degrees(LAM/RV):.2f}°  ({np.degrees(LAM/RV)*60:.0f}′)"],
 ["ring mode it maps to (m = 2π·14.3 Gpc / period)","—",f"m ≈ {circ/LAM:.0f}"],
 ["rms contrast σ at this scale (linear, σ₈=0.81)","0.19 (R = 37 h⁻¹ Mpc)","0.08 (R = 74 h⁻¹ Mpc)"]]
t=ax.table(cellText=rows,loc="center",cellLoc="center",colWidths=[0.5,0.25,0.25]); t.auto_set_font_size(False); t.set_fontsize(9.5); t.scale(1,1.9)
ax.set_title("(d) How the cell appears at the ring distance (ring sampled at 1°, i.e. m ≤ 180, so the cell itself is not resolved)",fontsize=10)
fig.text(0.5,0.06,"Period from Einasto et al. 1997 (arXiv:astro-ph/9701018): λ₀ = 120 ± 15 h⁻¹ Mpc → 171 Mpc for h = 0.7. Their picture: voids of almost equal size, rich superclusters in chains between them;\npoor superclusters and isolated clusters in the void walls. Void radius in (a) is our illustrative choice.",ha="center",fontsize=8)
fig.savefig("fig2_cell_and_scaling.png",dpi=130,bbox_inches="tight")
