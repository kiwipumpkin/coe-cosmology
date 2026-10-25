import sys,json,numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
ME='/workspace/coe-cosmology-gh/draft/current/analyses/void_mass_network';sys.path.insert(0,ME);import vmn_gen as G
g=G.generate(Mmin=3e15,seed=2000,comp_mode="lcdm");S=g["S"]/1000;XC=g["XC"]/1000;Rv=g["Rv"]/1000;cen=g["cen"]/1000
vv=np.linalg.norm(g["S"],axis=1)<g["params"]["R_out"];nv=int(vv.sum());nc=len(XC)
K=np.load('thin4_keep.npz');ku=K['ku'];r=json.load(open('thin4_results.json'))

FR=['0.01','0.0125','0.015','0.0175','0.02'];info={}
def sel(f):
    L=r['thin'][f];kb=int(np.argmin([d['chi2'] for d in L]));c2=L[kb]['chi2'];ff=float(f)
    rng=np.random.default_rng(20261021)
    for k in range(kb+1):uv=rng.random(nv);uc=rng.random(nc)
    return kb,c2,ku[:nv]&(uv<ff),ku[nv:]&(uc<ff)
CC=[-0.041,-0.035,0.1353,0.0];Rb=1.5*(1+1/2.937);aC=np.array([1,0,np.tan(np.radians(15))]);aC/=np.linalg.norm(aC);PC=42*aC
lnr=lambda P:sum(c*(np.minimum(np.linalg.norm(P-PC,axis=-1)/42,Rb)-1)**(k+1) for k,c in enumerate(CC))
def draw(ax,lim,hz,small,vidx,kc):
    xs=np.linspace(-lim,lim,300);XX,YY=np.meshgrid(xs,xs);Z=lnr(np.stack([XX,YY,np.zeros_like(XX)],-1))
    im=ax.pcolormesh(XX,YY,100*(np.exp(Z)-1),cmap='Greys',shading='auto',zorder=0);cs=ax.contour(XX,YY,100*(np.exp(Z)-1),8,colors='orange',linewidths=0.7,zorder=0.5);ax.clabel(cs,fmt='%.2f%%',fontsize=7)
    sv=vidx[np.abs(S[vidx,2])<hz];rc=np.sqrt(np.maximum(Rv[sv]**2-S[sv,2]**2,0))
    ax.add_collection(PatchCollection([plt.Circle(S[i,:2],q) for i,q in zip(sv,rc) if q>0],facecolor="#2471a3",alpha=0.45,edgecolor="none",zorder=1))
    sc=kc&(np.abs(XC[:,2])<hz);m=g["MC"][sc]
    ax.scatter(XC[sc,0],XC[sc,1],s=(3 if not small else 10)+(30 if not small else 80)*(m/g["MC"].max())**1.5,c="#c0392b",edgecolors="none",zorder=3)
    ax.add_patch(plt.Circle((0,0),14.3,fill=False,color="k",lw=1.4,zorder=4));ax.plot(0,0,"k+",ms=12,mew=2,zorder=5)
    ax.add_patch(plt.Circle((0,0),0.5,fill=False,color="k",lw=0.8,ls=":",zorder=4))
    ax.add_patch(plt.Circle((0,0),4,fill=False,color="tab:purple",lw=2.2,ls="--",zorder=4,label="4 Gpc sampling sphere"))
    ax.annotate('',xy=(0.85*lim*aC[0]/np.hypot(aC[0],aC[1]),0),xytext=(0,0),arrowprops=dict(arrowstyle='->',color='darkgreen',lw=2.5),zorder=6)
    ax.text(0.55*lim,0.05*lim,'to CoE (l=0°, b=15°), 42 Gpc',color='darkgreen',fontsize=11,zorder=6)
    ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim);ax.set_aspect("equal");ax.set_xlabel("Gpc (galactic x)");ax.set_ylabel("Gpc (galactic y)");plt.colorbar(im,ax=ax,shrink=0.7,label='ρ/ρ(us) − 1  (%) at z=0 plane')
    return len(sv),int(sc.sum())

pf,pa=plt.subplots(2,3,figsize=(24,15));pa=pa.ravel()
for i,f in enumerate(FR):
    kb,c2,kv,kc=sel(f);vidx=np.where(vv)[0][kv]
    fig,ax=plt.subplots(1,2,figsize=(22,10))
    a,b=draw(ax[0],g["params"]["R_out"]/1000+0.5,0.4,False,vidx,kc);draw(ax[1],6,0.4,True,vidx,kc)
    ax[0].set_title(f"Density-shaped thinning, kept f={f}, best draw {kb}, CMB χ²={c2:.1f} (web read at 4 Gpc)\nslab |z|<0.4 Gpc: {a} voids (blue), {b} clusters (red); solid circle 14.3 Gpc, + = us",fontsize=11)
    ax[1].set_title("Zoom 12×12 Gpc; purple dashed = 4 Gpc sampling sphere, dotted = 0.5 Gpc empty zone",fontsize=11)
    fig.tight_layout();fig.savefig(f"thin4_bodies_f{f}.png",dpi=90);plt.close(fig)
    draw(pa[i],6,0.4,True,vidx,kc);pa[i].set_title(f"f={f}, draw {kb}, χ²={c2:.1f}, slab voids/clusters {a}/{b}",fontsize=11)
    info[f]=(kb,c2,a,b,int(kv.sum()),int(kc.sum()))
pa[5].axis('off');pf.suptitle("Best-of-5 density-shaped thinned webs (4 Gpc read), |z|<0.4 Gpc slab, 12×12 Gpc zoom; purple dashed = 4 Gpc sampling sphere",fontsize=13)
pf.tight_layout();pf.savefig("thin4_bodies_panel.png",dpi=80);print(info)
