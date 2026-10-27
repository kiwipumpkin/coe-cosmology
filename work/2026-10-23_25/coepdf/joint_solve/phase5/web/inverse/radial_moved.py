# Radial distribution (r/R) of displaced bodies: original vs final positions, against the uniform-in-volume expectation; plus displacement histogram.
import glob,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
W='/workspace/coepdf/joint_solve/phase5/web/inverse/';R=21.7439;fs=sorted(glob.glob(W+'inv_s*.npz'))
fig,ax=plt.subplots(1,2,figsize=(13,5));e=np.linspace(0,1,21)
for f in fs:
    z=np.load(f);m=z['moved'];s=f.split('_')[-1][:-4]
    for a in ('y0','y','Resid'):assert np.all(np.isfinite(z[a])),a
    ax[0].hist(np.linalg.norm(z['y'][m],axis=1)/R,e,histtype='step',lw=1.5,label=f'{s} final (n={m.sum()})')
    ax[0].hist(np.linalg.norm(z['y0'][m],axis=1)/R,e,histtype='step',ls='--',label=f'{s} original')
    ax[1].hist(np.linalg.norm(z['y'][m]-z['y0'][m],axis=1),30,histtype='step',label=s)
x=np.linspace(0.01,0.975,200);n=int(np.mean([np.load(f)['moved'].sum() for f in fs]))
ax[0].plot(x,n*3*x**2/(0.975**3-0.01**3)*(e[1]-e[0]),'k:',label='uniform in volume')
ax[0].axvline(0.975,color='gray',lw=0.8);ax[0].set_xlabel('r / R (R = 21.74 Gpc shell)');ax[0].set_ylabel('moved bodies per bin');ax[0].legend(fontsize=7)
ax[1].set_xlabel('displacement [Gpc]');ax[1].set_ylabel('count');ax[1].legend(fontsize=7)
fig.suptitle('Inverse problem: where the minimally-displaced bodies come from / go');fig.tight_layout();fig.savefig(W+'radial_moved.png',dpi=130)
