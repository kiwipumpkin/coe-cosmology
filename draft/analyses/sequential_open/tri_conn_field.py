"""foam s2, full Delaunay triangulation of the slab: (1) degree + link-length histograms; (2) continuous density-contrast field.
Field on the slab mid-plane z=0, 0.05 Gpc grid: clusters = 3D Gaussian of width sigma = body radius (0.158 Gpc), delta = M/(rho_local (2 pi)^1.5 sigma^3) e^{-d^2/2sigma^2};
voids = the model's top-hat deficit delta_V = -0.70 (seq_core.DELTA_V) inside R, softened with a tanh edge of width 0.1 R. d = 3D distance from the body to the plane point."""
import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.spatial import Delaunay
src=open("tri_view.py").read(); exec(src[:src.index("js,pre=sys.argv")])
import seq_core
d=json.load(open("minedit_ent/b_foam_rim_s2_tidalcore_entropy.json")); base=d["tag"].replace("_entropy","")
LAY={"before":EB.layout(base),"after":EB.layout(base,d)}; CB="0.55"; CA="#1f3a93"; info={}
# ---------- (1) connections
T={}
for k,(X,r,v,M,_) in LAY.items():
    sl,P,S,E,A=tri(X,r); deg=np.bincount(E.ravel(),minlength=len(sl)); X3=X[sl]; L=np.linalg.norm(X3[E[:,0]]-X3[E[:,1]],axis=1); vs=v[sl]
    ty=np.where(vs[E[:,0]]&vs[E[:,1]],"vv",np.where(~vs[E[:,0]]&~vs[E[:,1]],"cc","cv")); T[k]=(deg,L,ty)
    info[k]=dict(nodes=len(sl),links=len(E),mean_degree=float(deg.mean()),median_link_Gpc=float(np.median(L)),**{f"median_{t}":float(np.median(L[ty==t])) for t in ("cc","cv","vv")},**{f"n_{t}":int((ty==t).sum()) for t in ("cc","cv","vv")})
fig,ax=plt.subplots(3,1,figsize=(8.5,13)); fig.subplots_adjust(hspace=0.62)
bi=np.arange(0,max(T["before"][0].max(),T["after"][0].max())+2)-0.5
for k,c in (("before",CB),("after",CA)):
    ax[0].hist(T[k][0],bi,histtype="step",color=c,lw=2.2,label=f"{k.upper()}: mean {T[k][0].mean():.2f} connections per node ({len(T[k][0]):,} nodes)")
    ax[1].hist(T[k][1],np.linspace(0,8,49),histtype="step",color=c,lw=2.2,label=f"{k.upper()}: median {np.median(T[k][1]):.2f} Gpc ({len(T[k][1]):,} links)")
ax[0].set_xlabel("connections per node (degree)",fontsize=12); ax[0].set_title("(a) Connections per node",fontsize=13)
ax[1].set_xlabel("link length, 3D separation (Gpc; links > 8 Gpc not shown)",fontsize=12); ax[1].set_xlim(0,8); ax[1].set_title("(b) Link lengths",fontsize=13)
x=np.arange(3)
for o,k,c in ((-0.2,"before",CB),(0.2,"after",CA)):
    vals=[info[k][f"median_{t}"] for t in ("cc","cv","vv")]; ax[2].bar(x+o,vals,0.38,color=c,label=k.upper())
    for xi,vv,t in zip(x,vals,("cc","cv","vv")): ax[2].text(xi+o,vv,f"{vv:.2f}\nn={info[k]['n_'+t]}",ha="center",va="bottom",fontsize=9.5)
ax[2].set_xticks(x); ax[2].set_xticklabels(["cluster–cluster","cluster–void","void–void"],fontsize=11); ax[2].set_ylabel("median link length (Gpc)",fontsize=12)
ax[2].set_title("(c) Median link length by link type",fontsize=13); ax[2].set_ylim(0,max(info[k][f"median_{t}"] for k in info for t in ("cc","cv","vv"))*1.35)
for a_ in ax[:2]: a_.set_ylabel("count",fontsize=12)
for a_ in ax: a_.tick_params(labelsize=11); a_.grid(alpha=0.3); a_.legend(loc="upper center",bbox_to_anchor=(0.5,-0.2),ncol=1 if a_ is not ax[2] else 2,frameon=False,fontsize=11)
fig.suptitle("foam s2, full Delaunay triangulation of the slab: BEFORE vs AFTER",fontsize=13.5,fontweight="bold",y=0.93)
fig.savefig("tri_foam_s2_conn_hist.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# ---------- (2) field
h=0.05; g=np.arange(-30,30+h/2,h); GX,GY=np.meshgrid(g,g); disc=np.hypot(GX,GY)<=30
def field(X,r,v,M):
    rho=seq_core.rho_gpc(X)*seq_core.GPC**3; F=np.zeros(GX.shape)   # kg per Gpc^3
    for i in np.where((r>0)&(np.linalg.norm(X,axis=1)<=31))[0]:
        if v[i]:
            R=r[i]; ext=R*1.3
            if abs(X[i,2])>ext: continue
        else:
            R=r[i]; ext=4*R
            if abs(X[i,2])>ext: continue
        i0=np.searchsorted(g,X[i,0]-ext); i1=np.searchsorted(g,X[i,0]+ext); j0=np.searchsorted(g,X[i,1]-ext); j1=np.searchsorted(g,X[i,1]+ext)
        dx=GX[j0:j1,i0:i1]-X[i,0]; dy=GY[j0:j1,i0:i1]-X[i,1]; dd=np.sqrt(dx**2+dy**2+X[i,2]**2)
        if v[i]: F[j0:j1,i0:i1]+=seq_core.DELTA_V*0.5*(1-np.tanh((dd-R)/(0.1*R)))
        else: F[j0:j1,i0:i1]+=M[i]/(rho[i]*(2*np.pi)**1.5*R**3)*np.exp(-dd**2/(2*R**2))
    F[~disc]=np.nan; return F
FLD={k:field(L[0],L[1],L[2],L[3]) for k,L in LAY.items()}
ring=disc&(np.abs(np.hypot(GX,GY)-RMIN)<3)
for k,F in FLD.items():
    comp=np.sign(F)*np.log10(1+np.abs(F))
    info[k].update(field_rms_all=float(np.sqrt(np.nanmean(F**2))),field_rms_ring3=float(np.sqrt(np.nanmean(F[ring]**2))),field_mean_all=float(np.nanmean(F)),field_mean_ring3=float(np.nanmean(F[ring])),
                   field_max=float(np.nanmax(F)),field_min=float(np.nanmin(F)),void_frac_all=float(np.nanmean(F[disc]<-0.3)),void_frac_ring3=float(np.mean(F[ring]<-0.3)),
                   comp_rms_all=float(np.sqrt(np.nanmean(comp**2))),comp_rms_ring3=float(np.sqrt(np.nanmean(comp[ring]**2))))
json.dump(info,open("tri_foam_s2_conn_field.json","w"),indent=1); print(json.dumps(info,indent=1))
HT=lambda F: np.sign(F)*np.log10(1+np.abs(F))           # display height: signed log10(1+|delta|)
zmax=max(np.nanmax(HT(F)) for F in FLD.values()); zmin=min(np.nanmin(HT(F)) for F in FLD.values()); norm=TwoSlopeNorm(vmin=-max(abs(zmin),0.3)*1.0,vcenter=0,vmax=zmax)
cmap=plt.get_cmap("RdBu_r")
# (a) surface: field smoothed with a 0.5 Gpc Gaussian for display (cluster bumps are 0.16 Gpc wide = needles at this scale), linear delta
from scipy.ndimage import gaussian_filter
SM={}
for k,F in FLD.items():
    F0=np.nan_to_num(F,nan=0.0); SM[k]=gaussian_filter(F0,0.5/h); SM[k][~disc]=np.nan
    info[k].update(smooth05_rms_all=float(np.sqrt(np.nanmean(SM[k]**2))),smooth05_rms_ring3=float(np.sqrt(np.nanmean(SM[k][ring]**2))),smooth05_max=float(np.nanmax(SM[k])),smooth05_min=float(np.nanmin(SM[k])))
json.dump(info,open("tri_foam_s2_conn_field.json","w"),indent=1)
ZC=2.0; nS=TwoSlopeNorm(vmin=-0.8,vcenter=0,vmax=ZC); lim=ZC   # peaks clipped at delta = 2 for display
fig=plt.figure(figsize=(15,7.8))
for n,k in enumerate(("before","after")):
    a_=fig.add_subplot(1,2,n+1,projection="3d"); s_=4; Z=np.clip(np.nan_to_num(SM[k][::s_,::s_],nan=0.0),-0.8,ZC); Xs=GX[::s_,::s_]; Ys=GY[::s_,::s_]
    a_.plot_surface(Xs,Ys,Z,facecolors=cmap(nS(Z)),rstride=1,cstride=1,linewidth=0,antialiased=False,shade=False)
    th=np.linspace(0,2*np.pi,300); a_.plot(RMIN*np.cos(th),RMIN*np.sin(th),0*th+ZC,color="k",lw=1.6,zorder=10)
    u=C_G[:2]/np.linalg.norm(C_G[:2]); a_.quiver(0,0,ZC,30*u[0],30*u[1],0,color="#a93226",lw=2,arrow_length_ratio=0.08); a_.text(24*u[0],24*u[1],ZC*1.05,"to CoE",color="#a93226",fontsize=10,fontweight="bold")
    a_.set_zlim(-0.8,ZC); a_.view_init(elev=35,azim=-60); a_.set_xlabel("x (Gpc)",fontsize=11); a_.set_ylabel("y (Gpc)",fontsize=11); a_.set_zlabel("δ (0.5 Gpc smoothing)",fontsize=11)
    a_.set_title(f"{k.upper()}: rms δ {info[k]['smooth05_rms_all']:.3f} (all), {info[k]['smooth05_rms_ring3']:.3f} (≤3 Gpc of ring)",fontsize=12.5); a_.tick_params(labelsize=9)
sm=plt.cm.ScalarMappable(norm=nS,cmap=cmap); cb=fig.colorbar(sm,ax=fig.axes,orientation="horizontal",fraction=0.04,pad=0.04,aspect=60)
cb.set_label("density contrast δ, smoothed 0.5 Gpc, peaks clipped at δ = 2 (max ≈ 4.1)\nred = cluster peaks, blue = void troughs, white ≈ mean density; black circle = CMB ring and arrow = from us toward the CoE (both drawn on the top plane)",fontsize=11)
fig.suptitle("foam s2 — continuous density-contrast field in the slab mid-plane (same z-scale and view)",fontsize=13.5,fontweight="bold")
fig.savefig("tri_foam_s2_field_surface.png",dpi=130,bbox_inches="tight"); plt.close(fig)
# (b) flat map + residual
fig=plt.figure(figsize=(13,11.2)); GS=fig.add_gridspec(3,2,height_ratios=[1.0,0.035,0.42],hspace=0.42,wspace=0.1,left=0.06,right=0.98,top=0.9,bottom=0.07)
for n,k in enumerate(("before","after")):
    a_=fig.add_subplot(GS[0,n]); im=a_.pcolormesh(GX,GY,HT(FLD[k]),cmap=cmap,norm=norm,shading="auto",rasterized=True,zorder=0)
    a_.add_patch(Circle((0,0),RMIN,fill=False,ec="k",lw=2.0,zorder=6)); a_.plot(0,0,"k+",ms=11,mew=2,zorder=9); frame30(a_,lw=1.2,fs=EB.FK)
    a_.set_xticks([-30,0,30]); a_.set_yticks([-30,0,30]); a_.tick_params(labelsize=EB.FK)
    a_.set_title(f"{k.upper()}: rms δ {info[k]['field_rms_all']:.2f} overall, {info[k]['field_rms_ring3']:.2f} within 3 Gpc of ring",fontsize=EB.FT)
cb=fig.colorbar(im,cax=fig.add_subplot(GS[1,:]),orientation="horizontal"); cb.set_label("density contrast as sign·log₁₀(1+|δ|): red = clusters, blue = voids (δ = −0.70 → −0.23), white ≈ mean",fontsize=EB.FK); cb.ax.tick_params(labelsize=EB.FK)
EB.residual_below(fig,GS[2,:],d,LAY["before"][4])
fig.legend(handles=[Line2D([],[],color="k",lw=2,label="CMB ring / our sphere"),Line2D([],[],color="0.4",ls=":",lw=1.2,label="30 Gpc edge of model"),
    Line2D([],[],color="#a93226",lw=2,marker=">",label="vector from us toward the CoE"),Rectangle((0,0),1,1,fc=EB.BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)")],
    loc="upper center",ncol=4,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.03))
fig.suptitle("foam s2 — density-contrast field, slab mid-plane, BEFORE | AFTER",fontsize=EB.FT+1,fontweight="bold",color="#8e44ad")
fig.savefig("tri_foam_s2_field_map.png",dpi=140,bbox_inches="tight")
