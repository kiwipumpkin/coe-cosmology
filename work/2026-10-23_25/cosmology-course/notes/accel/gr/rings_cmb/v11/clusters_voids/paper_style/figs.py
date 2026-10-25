import os,sys,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.patches import Circle
SEQ='/workspace/clump_scan/sequential';c=os.getcwd();os.chdir(SEQ);sys.path.insert(0,SEQ);import seq_core as S;os.chdir(c)
f=np.load('/workspace/clump_scan/sequential_open/final_comp_s3.npz');X,r,M,g=f['Xp'],f['rp'],f['Mp'],f['grp'];V=f['kp']=='void';Cl=~V
nV,nC=int(V.sum()),int(Cl.sum());rho=S.rho_gpc(X);GPC=3.0857e25
dens=np.abs(M)/(4/3*np.pi*(r*GPC)**3)/rho if np.nanmax(rho)<1e-10 else np.abs(M)/(4/3*np.pi*r**3)/rho
CG=np.array(S.C_G);A=float(S.A_IN_G);print('C_G',CG,np.linalg.norm(CG),'core radius',A,'dens sample',dens[V][:3],dens[Cl][:3])
cap=f"Chris's compensated random web generator (compweb_gen.py), seed 3: {nV+nC:,} bodies = {nV:,} voids + {nC:,} clusters. Each void plus its wall clusters has zero net mass."
# coordinates: rotate so CoE is along +x in the slab plane
u=CG/np.linalg.norm(CG);w=np.cross(u,[0,0,1.]);w/=np.linalg.norm(w);z=np.cross(u,w)
P=np.c_[-(X@u),X@w,X@z]   # mirrored: CoE on the LEFT
# 1 side-view slab
sl=np.abs(P[:,2])<3.0
fig,ax=plt.subplots(figsize=(15,15))
for i in np.where(sl&V)[0]:ax.add_patch(Circle(P[i,:2],np.sqrt(max(r[i]**2-P[i,2]**2,0)),fc='#9ecae1',ec='#3182bd',alpha=.35,lw=.6))
ax.scatter(P[sl&V,0],P[sl&V,1],s=14,c='#08519c',label='void centre (light-blue disc = void cut by the slab)',zorder=3)
ms=M[sl&Cl]/np.median(M[Cl]);ax.scatter(P[sl&Cl,0],P[sl&Cl,1],s=10*ms,c='#cb181d',label='cluster (dot size ∝ mass)',zorder=4)
ax.add_patch(Circle((0,0),14.3,fill=False,ec='k',lw=2.5));ax.plot([],[],'k-',lw=2.5,label='our sphere / CMB ring, 14.3 Gpc')
ax.add_patch(Circle((-np.linalg.norm(CG),0),A,fc='pink',ec='m',alpha=.6,lw=1.5));ax.fill([],[],fc='pink',ec='m',label=f'CoE core cap (radius {A:.1f} Gpc, centre {np.linalg.norm(CG):.1f} Gpc away)')
ax.annotate('',xy=(-30,0),xytext=(0,0),arrowprops=dict(arrowstyle='-|>',color='#a93226',lw=3,mutation_scale=25));ax.text(-30,1.2,'to CoE (58.2 Gpc)',color='#a93226',fontsize=15,fontweight='bold')
ax.plot(0,0,'k*',ms=18,label='us');ax.add_patch(Circle((0,0),30,fill=False,ec='.5',ls=':',lw=1.2))
ax.set_xlim(-32,32);ax.set_ylim(-32,32);ax.set_aspect('equal');ax.set_xlabel('distance along the CoE direction (Gpc; CoE to the left)',fontsize=16);ax.set_ylabel('distance across (Gpc)',fontsize=16);ax.tick_params(labelsize=13)
ax.set_title('Random compensated web, one seed: side view, slab |height| < 3 Gpc',fontsize=18,pad=95)
ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),ncol=2,fontsize=13,frameon=False)
fig.text(.5,.05,cap+f'\nShown: {int((sl&V).sum())} voids and {int((sl&Cl).sum())} clusters with centres inside the slab. Dotted circle: 30 Gpc modelled region.',ha='center',fontsize=12)
fig.savefig('web_sideview_seed3.png',dpi=80,bbox_inches='tight');plt.close(fig)
# 2 zoom on a few voids, each group drawn in its own frame (projected around its centre)
gv=np.where(V)[0];Rv=r[gv];pick=gv[np.argsort(np.abs(Rv-np.median(Rv)))[:4]]
fig,axs=plt.subplots(2,2,figsize=(16,16))
for ax,iv in zip(axs.flat,pick):
    gi=g[iv];mem=np.where((g==gi)&Cl)[0];R=r[iv];d=P[mem]-P[iv];dist=np.linalg.norm(d,axis=1)
    ax.add_patch(Circle((0,0),R,fc='#9ecae1',ec='#3182bd',alpha=.4,lw=2));ax.add_patch(Circle((0,0),1.3*R,fill=False,ec='k',ls='--',lw=1.5))
    for j,dd in zip(mem,d):ax.add_patch(Circle(dd[:2],r[j],fc='#cb181d',ec='k',alpha=.8,lw=.5))
    # neighbours (other bodies) faintly
    nb=np.where((np.linalg.norm(X-X[iv],axis=1)<1.8*R)&(g!=gi))[0]
    for j in nb:ax.add_patch(Circle((P[j]-P[iv])[:2],r[j],fc='.85' if Cl[j] else 'none',ec='.6',lw=.5))
    Mdef=-M[iv];held=M[mem].sum();dv=-dens[iv]
    ax.set_title(f'void radius {R:.2f} Gpc, depth {dv:+.2f}\nmissing mass {Mdef:.2e} kg = {len(mem)} clusters holding {held:.2e} kg\nnet mass of void+wall {(held-Mdef)/Mdef:+.1e} of deficit; clusters at {dist.min()/R:.2f}–{dist.max()/R:.2f} void radii',fontsize=13)
    ax.set_xlim(-1.9*R,1.9*R);ax.set_ylim(-1.9*R,1.9*R);ax.set_aspect('equal');ax.tick_params(labelsize=12);ax.set_xlabel('Gpc (projected, from void centre)',fontsize=13)
h=[plt.Line2D([],[],marker='o',ls='',mfc='#9ecae1',mec='#3182bd',ms=16),plt.Line2D([],[],ls='--',c='k'),plt.Line2D([],[],marker='o',ls='',mfc='#cb181d',mec='k',ms=10),plt.Line2D([],[],marker='o',ls='',mfc='.85',mec='.6',ms=10)]
fig.legend(h,['void (light blue disc = its radius)','wall edge at 1.3 void radii','its own compensating clusters (true size)','other bodies nearby (grey)'],loc='upper center',ncol=2,fontsize=14,frameon=False,bbox_to_anchor=(.5,1.03))
fig.suptitle('Compensation rule: each void’s missing mass sits in clusters in its wall (1–1.3 void radii)',fontsize=18,y=1.07)
fig.text(.5,-.01,cap+' Four voids of median size; clusters projected onto the page, so some appear inside the disc.',ha='center',fontsize=12)
fig.tight_layout();fig.savefig('compensation_zoom_seed3.png',dpi=80,bbox_inches='tight');plt.close(fig)
# 3 distributions
fig,axs=plt.subplots(1,3,figsize=(22,7.5))
a=axs[0];a.hist(r[V],bins=50,color='#3182bd');q=np.percentile(r[V],[5,50,95]);a.set_title(f'void radii: median {q[1]:.2f} Gpc, 5–95% {q[0]:.2f}–{q[2]:.2f} Gpc',fontsize=14);a.set_xlabel('void radius (Gpc)',fontsize=14)
a=axs[1];a.hist(M[Cl],bins=np.geomspace(M[Cl].min(),M[Cl].max(),50),color='#cb181d');a.set_xscale('log');a.set_title(f'cluster masses: median {np.median(M[Cl]):.2e} kg',fontsize=14);a.set_xlabel('cluster mass (kg)',fontsize=14)
a=axs[2];a.hist(dens[Cl],bins=np.geomspace(25,350,40),color='#cb181d',alpha=.8,label=f'cluster contrast (median {np.median(dens[Cl]):.0f}; median radius {np.median(r[Cl]):.3f} Gpc)');a.set_xscale('log');a.set_xlabel('density contrast',fontsize=14)
a2=a.twiny();a2.hist(-dens[V],bins=30,color='#3182bd',alpha=.6,label='void depth (top axis)');a2.set_xlabel('void depth (fractional density contrast)',fontsize=13)
a.set_title('cluster contrast (bottom) and void depth (top)',fontsize=14,pad=45);a.legend(loc='upper left',fontsize=11);a2.legend(loc='upper right',fontsize=11)
for a in axs:a.set_ylabel('number',fontsize=14);a.tick_params(labelsize=12)
fig.suptitle('Input laws as drawn in this seed: Jennings–Li–Hu void sizes, Tinker cluster masses, uniform void depths, log-uniform cluster contrasts',fontsize=16,y=1.05)
fig.text(.5,-.04,cap+f' Void volume fraction target 29–30%.',ha='center',fontsize=12)
fig.tight_layout();fig.savefig('input_distributions_seed3.png',dpi=80,bbox_inches='tight')
print(nV,nC,q,np.median(M[Cl]),np.median(dens[Cl]),np.median(r[Cl]),np.percentile(-dens[V],[0,100]))
