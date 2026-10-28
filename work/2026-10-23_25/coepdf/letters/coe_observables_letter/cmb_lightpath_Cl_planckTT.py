# D_l of the light-path total map (monopole removed; and monopole+dipole removed) vs Planck TT R3.01, l<=767 (3*Nside-1).
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
m=np.load('cmb_lightpath_maps.npz');T=2.7255e6;p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt')
l=np.arange(len(m['cl']));f=l*(l+1)/2/np.pi*T**2
fig,ax=plt.subplots(figsize=(9,5.5))
ax.loglog(p[:,0],p[:,1],'k',lw=1,label='Planck TT (R3.01)')
ax.loglog(l[1:],f[1:]*m['cl'][1:],'C3',label='model total, monopole removed')
ax.loglog(l[2:],f[2:]*m['cl_nodip'][2:],'C0',label='model total, monopole+dipole removed')
ax.set_xlabel('$\\ell$');ax.set_ylabel('$D_\\ell$ [$\\mu$K$^2$]');ax.legend();ax.set_title('Light-path prediction (R_see run, Nside 256), unscaled')
plt.savefig('cmb_lightpath_Cl_planckTT.png',dpi=120,bbox_inches='tight')
