# v3: radial distribution r/R of CHANGED bodies (original vs final) vs uniform-in-volume, v2 seed-1 moved-final for comparison;
# displacement histogram (position-moved bodies) and delta change histogram. Usage: python3 radial.py [dir/]. Fails on NaN.
import sys,glob,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
W='/workspace/coepdf/joint_solve/phase5/web/inverse/v3/'+(sys.argv[1] if len(sys.argv)>1 else '');R=21.7439;fs=sorted(glob.glob(W+'inv3_s*.npz'))
fig,ax=plt.subplots(1,3,figsize=(17,5));e=np.linspace(0,1,21)
for f in fs:
    z=np.load(f);c=z['changed'];mv=z['moved'];s=f.split('_')[-1][:-4]
    for k in ('y0','y','delta','Resid'):
        if not np.all(np.isfinite(z[k])):raise SystemExit('NONFINITE '+k)
    ax[0].hist(np.linalg.norm(z['y'][c],axis=1)/R,e,histtype='step',lw=1.5,label=f'v3 {s} final (n={c.sum()})');ax[0].hist(np.linalg.norm(z['y0'][c],axis=1)/R,e,histtype='step',ls='--',label=f'v3 {s} original')
    ax[1].hist(np.linalg.norm(z['y'][mv]-z['y0'][mv],axis=1),30,histtype='step',label=f'v3 {s} (n={mv.sum()})');ax[2].hist((z['delta']-z['delta0'])[c],40,histtype='step',label=f'v3 {s}')
v=np.load('/workspace/coepdf/joint_solve/phase5/web/inverse/inv_s1.npz');ax[0].hist(np.linalg.norm(v['y'][v['moved']],axis=1)/R,e,histtype='step',color='0.5',lw=1,label='v2 s1 final (positions only)')
x=np.linspace(0.01,0.975,200);n=int(np.mean([np.load(f)['changed'].sum() for f in fs]))
ax[0].plot(x,n*3*x**2/(0.975**3-0.01**3)*(e[1]-e[0]),'k:',label='uniform in volume');ax[0].axvline(0.975,color='gray',lw=0.8)
ax[0].set_xlabel('r / R (R = 21.74 Gpc shell; line = 0.975 R geometric limit)');ax[0].set_ylabel('changed bodies per bin');ax[0].legend(fontsize=7)
ax[1].set_xlabel('displacement [Gpc]');ax[1].set_ylabel('count');ax[1].legend(fontsize=7);ax[2].set_yscale('log')
ax[2].set_xlabel('delta_final - delta_start (changed bodies)');ax[2].set_ylabel('count');ax[2].legend(fontsize=7)
fig.suptitle('v3 inverse web: where changed bodies come from / go, how far, how much mass changes');fig.tight_layout();fig.savefig(W+'radial.png',dpi=130)
