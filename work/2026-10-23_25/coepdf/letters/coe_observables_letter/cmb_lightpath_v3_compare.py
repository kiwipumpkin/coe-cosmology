# v3 comparison: monopole+dipole-removed maps for rules A, B, C (uK, each own colour scale) and their D_l vs Planck TT.
import numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from cmb_lightpath_v3_common import *
fig=plt.figure(figsize=(16,9));L,Dl=planck();ax=fig.add_subplot(2,1,2)
for i,R in enumerate('ABC'):
    m,m1,cl,cl1,kind,mean=stats(R);hp.mollview(m1*T0uK,fig=fig.number,sub=(2,3,i+1),title=f'v3{R} mono+dipole removed [uK], rms {np.std(m1)*T0uK:.3g}',cmap='coolwarm')
    l=np.arange(len(cl1));ax.loglog(l[2:],l[2:]*(l[2:]+1)/2/np.pi*cl1[2:]*T0uK**2,label=f'v3{R}')
ax.loglog(L,Dl,'k',label='Planck TT');ax.set_xlabel('l');ax.set_ylabel('D_l [uK^2] (mono+dipole removed)');ax.legend()
plt.savefig('cmb_lightpath_v3_compare.png',dpi=100,bbox_inches='tight')
