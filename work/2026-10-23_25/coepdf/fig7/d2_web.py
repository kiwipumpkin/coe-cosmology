# REGENERATED for rev7 (replaces rev5/6 Press-Schechter/void-function web ingredients, which used LCDM sigma8/Omega_m inputs):
# the compensating bodies actually used by the solve: p5_c_p8 (final run: background D=34.96 plus (1-f)=0.74 of the p8 field l=2 on the seed-42 sky, bodies inside R=R_ls=0.18957 Gpc, |M|<=1e14 Msun each).
import numpy as np,json,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
P='/workspace/coepdf/joint_solve/phase2/';f='p5_c_p8';Z=np.load(P+f+'.npz');J=json.load(open(P+f+'.json'));R=J['R']
y,M=Z['y'],Z['M'];r=np.linalg.norm(y,axis=1);c=M>0;v=M<0
fig=plt.figure(figsize=(13,4.3))
a=fig.add_subplot(1,3,1);bins=np.geomspace(1e10,1e14,25);a.hist(M[c],bins,color='b',alpha=.6,label=f'clusters (M>0): {c.sum()}');a.hist(-M[v],bins,color='r',alpha=.6,label=f'voids (M<0): {v.sum()}')
a.set_xscale('log');a.set_xlabel('|M| [M$_\\odot$]');a.set_ylabel('count');a.legend(fontsize=8);a.set_title(f"$\\Sigma M_+$={M[c].sum():.2e}, $\\Sigma M_-$={M[v].sum():.2e}",fontsize=9)
a=fig.add_subplot(1,3,2);a.hist(r[c]/R,30,color='b',alpha=.6,weights=M[c],label='clusters (mass-weighted)');a.hist(r[v]/R,30,color='r',alpha=.6,weights=-M[v],label='voids (|mass|-weighted)')
a.set_xlabel('distance from observer / R');a.set_ylabel('M$_\\odot$');a.legend(fontsize=8)
hp.mollview(Z['C']*1e6,sub=(1,3,3),title='bodies-only sky $\\Delta T/T\\times10^6$',fig=fig.number,cbar=True)
fig.tight_layout();fig.savefig('d2_web.png',dpi=130,facecolor='white')
