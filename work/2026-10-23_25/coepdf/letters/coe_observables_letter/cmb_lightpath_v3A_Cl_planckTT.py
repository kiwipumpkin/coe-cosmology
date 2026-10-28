# v3A: D_l of total (monopole removed; monopole+dipole removed) vs Planck TT R3.01, l<=767, unscaled.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from cmb_lightpath_v3_common import *
m,m1,cl,cl1,kind,mean=stats('A');l=np.arange(len(cl));f=l*(l+1)/2/np.pi*T0uK**2;L,Dl=planck()
fig,ax=plt.subplots(figsize=(9,5.5));ax.loglog(L,Dl,'k',lw=1,label='Planck TT');ax.loglog(l[1:],f[1:]*cl[1:],'C3',label='monopole removed');ax.loglog(l[2:],f[2:]*cl1[2:],'C0',label='monopole+dipole removed')
ax.set_xlabel('l');ax.set_ylabel('D_l [uK^2]');ax.legend();ax.set_title('v3A light-path prediction, Nside 256')
plt.savefig('cmb_lightpath_v3A_Cl_planckTT.png',dpi=110,bbox_inches='tight')
