# rev7 regeneration of rev6 Fig 11 (cluster-void web). rev6 used vmn_gen seed 2000 in the D=42 era domain; rev7 draws the actual
# 3000 bodies of the converged CMB solve (joint_solve/phase2/p5_c_p8.npz: positions y [Gpc], signed masses M) inside R_ls.
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
z=np.load('/workspace/coepdf/joint_solve/phase2/p5_c_p8.npz');y,M=z['y'],z['M'];R=0.18957;r=np.linalg.norm(y,axis=1)
ins=r<=R;cl=ins&(M>0);vo=ins&(M<0);hw=0.02
fig,ax=plt.subplots(1,3,figsize=(20,6.6))
for a,(lim,ttl) in zip(ax[:2],((R*1.05,f'slab |z|<{hw} Gpc, full sphere R_ls={R} Gpc'),(R*0.35,'zoom (inner 0.066 Gpc)'))):
    s=np.abs(y[:,2])<hw;sz=lambda m:6+60*(np.abs(m)/np.abs(M).max())**0.5
    a.scatter(*y[s&cl][:,:2].T,s=sz(M[s&cl]),c='crimson',alpha=.7,label=f'clusters (M>0): {int((s&cl).sum())} in slab')
    a.scatter(*y[s&vo][:,:2].T,s=sz(M[s&vo]),c='royalblue',alpha=.7,label=f'voids (M<0): {int((s&vo).sum())} in slab')
    a.add_patch(plt.Circle((0,0),R,fill=False,ls='--',color='k'));a.set_xlim(-lim,lim);a.set_ylim(-lim,lim);a.set_aspect('equal')
    a.set_xlabel('x [Gpc]');a.set_ylabel('y [Gpc]');a.set_title(ttl,fontsize=10);a.legend(fontsize=8,loc='upper right')
b=np.linspace(0,R,25);rc=(b[1:]+b[:-1])/2;V=4/3*np.pi*(b[1:]**3-b[:-1]**3)
nc=np.histogram(r[cl],b)[0]/V;nv=np.histogram(r[vo],b)[0]/V;mc=np.histogram(r[cl],b,weights=M[cl])[0];mv=np.histogram(r[vo],b,weights=M[vo])[0]
ax[2].plot(rc,nc,'r-o',ms=3,label='cluster number density');ax[2].plot(rc,nv,'b-o',ms=3,label='void number density');ax[2].set_yscale('log')
ax[2].set_xlabel('radius from observer [Gpc]');ax[2].set_ylabel('bodies / Gpc³');ax[2].set_title(f'{int(cl.sum())} clusters, {int(vo.sum())} voids inside R_ls; net mass {M[ins].sum():.2e} M☉',fontsize=10);ax[2].legend(fontsize=8)
fig.tight_layout();fig.savefig(__file__.replace('.py','.png'),dpi=95);json.dump(dict(n_cl=int(cl.sum()),n_vo=int(vo.sum()),Mnet=float(M[ins].sum()),Mpos=float(M[cl].sum()),Mneg=float(M[vo].sum())),open(__file__.replace('.py','.json'),'w'));print(cl.sum(),vo.sum(),M[ins].sum())
