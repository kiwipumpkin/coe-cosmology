# rev7 regeneration of rev6 Fig 12 (smoothing ladder of the web density). rev6 smoothed a D=42-era vmn_gen field on Gpc kernels; rev7
# deposits the solved bodies (p5_c_p8.npz) on a 128^3 grid over the R_ls sphere and Gaussian-smooths with kernels scaled to R_ls.
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
z=np.load('/workspace/coepdf/joint_solve/phase2/p5_c_p8.npz');y,M=z['y'],z['M'];R=0.18957;N=128;L=2*R;dx=L/N
H=np.histogramdd(y,bins=[np.linspace(-R,R,N+1)]*3,weights=M)[0]/dx**3;c=(np.arange(N)+.5)*dx-R;X,Y,Z=np.meshgrid(c,c,c,indexing='ij');ins=X**2+Y**2+Z**2<=R**2
ks=[0.005,0.01,0.02,0.04,0.08];fig,ax=plt.subplots(1,len(ks)+1,figsize=(27,5));out={}
for a,k in zip(ax,ks):
    F=gaussian_filter(H,k/dx,mode='constant');Fs=F[:,:,N//2].copy();Fs[~ins[:,:,N//2]]=np.nan;v=np.nanpercentile(np.abs(Fs),99)
    a.imshow(Fs.T,origin='lower',extent=[-R,R,-R,R],cmap='RdBu_r',vmin=-v,vmax=v);a.set_title(f'σ_s={k} Gpc ({k/R:.3f} R_ls)',fontsize=10);a.set_xlabel('x [Gpc]')
    out[k]=dict(rms=float(F[ins].std()),pos_frac=float((F[ins]>0).mean()))
ax[0].set_ylabel('y [Gpc]  (slab z=0)')
ax[-1].loglog(ks,[out[k]['rms'] for k in ks],'ko-');ax[-1].set_xlabel('smoothing σ_s [Gpc]');ax[-1].set_ylabel('rms signed mass density [M☉/Gpc³]');ax[-1].set_title('ladder: rms vs scale',fontsize=10);ax[-1].grid(alpha=.3)
fig.tight_layout();fig.savefig(__file__.replace('.py','.png'),dpi=85);json.dump(out,open(__file__.replace('.py','.json'),'w'),indent=1);print(out)
