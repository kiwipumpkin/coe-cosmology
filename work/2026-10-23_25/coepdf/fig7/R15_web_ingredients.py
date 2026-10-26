# rev7 script for rev6 Fig 15 (web ingredients), computed from book/code/26_cosmo.py (Eisenstein-Hu transfer, Press-Schechter halo
# mass function, Sheth & van de Weygaert void size function). These are standard-cosmology inputs, independent of the CoE background.
# Right: compensated void placement - a void of radius R_v with interior contrast delta_v=-0.8 and a ridge (R_v..1.5R_v) holding
# exactly the missing mass, so the enclosed mass excess returns to zero at 1.5 R_v (computed numerically).
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
import sys,importlib.util;sp=importlib.util.spec_from_file_location('cosmo','/workspace/coepdf/book/code/26_cosmo.py');C=importlib.util.module_from_spec(sp);sp.loader.exec_module(C);ps_dndlnM=lambda M:C.ps_dndlnM(M)[0];svdw_dndlnR=C.svdw_dndlnR
M=np.geomspace(1e12,3e15,80);R=np.geomspace(2,60,80);fig,ax=plt.subplots(1,3,figsize=(17,4.8))
ax[0].loglog(M,ps_dndlnM(M),'k');ax[0].set_xlabel('M [M☉]');ax[0].set_ylabel('dn/dlnM [Mpc⁻³]');ax[0].set_title('Press–Schechter halo mass function (EH98 P(k), σ8=0.8)',fontsize=9)
ax[1].loglog(R,svdw_dndlnR(R),'b');ax[1].set_xlabel('void radius R_nl [Mpc]');ax[1].set_ylabel('dn/dlnR [Mpc⁻³]');ax[1].set_title('Sheth & van de Weygaert void size function (δ_v=−2.7)',fontsize=9)
x=np.linspace(0,2,2001);dv=-0.8;rid=-dv*1/(1.5**3-1);d=np.where(x<1,dv,np.where(x<1.5,rid,0.0));Menc=np.cumsum(d*3*x**2)*(x[1]-x[0])
ax[2].plot(x,d,'b',label='δ(r)');ax[2].plot(x,Menc,'r--',label='enclosed mass excess / (4π/3 ρ̄ R_v³)');ax[2].axhline(0,color='k',lw=.6)
ax[2].set_xlabel('r / R_v');ax[2].legend(fontsize=8);ax[2].set_title(f'compensated void: ridge δ={rid:.3f}, M_enc(1.5R_v)={Menc[np.searchsorted(x,1.5)]:+.4f}',fontsize=9)
fig.tight_layout();fig.savefig(__file__.replace('.py','.png'),dpi=120)
