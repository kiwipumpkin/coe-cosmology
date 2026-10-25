# REMADE (not original): axisymmetric harmonics l=0..7 on a ring + power bar
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.special import eval_legendre
rng=np.random.default_rng(0);th=np.linspace(0,2*np.pi,400,endpoint=False)
fig=plt.figure(figsize=(9.7,13.7));fig.suptitle('Mass-density angular harmonics and their power spectra\nBlue: overdensity | Orange: underdensity',fontsize=13)
for L in range(8):
  r,c=divmod(L,2);ax=fig.add_axes([0.03+c*0.5,0.765-r*0.24,0.27,0.19]);ax.set_title(f'ℓ = {L}, m = 0')
  v=eval_legendre(L,np.cos(th));rad=1+0.02*rng.standard_normal(th.size)
  ax.scatter(rad*np.sin(th),rad*np.cos(th),s=6,c=np.where(v>=0,'royalblue','darkorange'));ax.set_aspect('equal');ax.axis('off')
  b=fig.add_axes([0.31+c*0.5,0.785-r*0.24,0.15,0.105]);p=np.zeros(8);p[L]=1;b.bar(range(8),p,color='tab:blue')
  b.set_xticks(range(8));b.set_yticks([0,1]);b.set_ylim(0,1.15);b.tick_params(labelsize=6);b.set_xlabel('Multipole L',fontsize=6)
fig.savefig('f11_remade.png',dpi=110,facecolor='white')
