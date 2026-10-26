# rev7 script for rev6 Fig 10 (axisymmetric density harmonics and their power). Background-independent geometry (no D or R enters).
# Each panel: a ring of 400 bodies whose density sign follows P_l(cos θ) with 2% radial scatter; the bar chart is the power
# actually measured from those bodies by Legendre projection (a_L = (2L+1)/2 ∫ s(θ) P_L sinθ dθ), not set by hand.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.special import eval_legendre
rng=np.random.default_rng(0);th=np.linspace(0,2*np.pi,400,endpoint=False);mu=np.cos(th);t=np.linspace(0,np.pi,2001)
fig=plt.figure(figsize=(9.7,13.7));fig.suptitle('Mass-density angular harmonics and their measured power\nBlue: overdensity | Orange: underdensity',fontsize=13)
for L in range(8):
    r,c=divmod(L,2);ax=fig.add_axes([0.03+c*0.5,0.765-r*0.24,0.27,0.19]);ax.set_title(f'ℓ = {L}, m = 0')
    v=eval_legendre(L,mu);rad=1+0.02*rng.standard_normal(th.size)
    ax.scatter(rad*np.sin(th),rad*np.cos(th),s=6,c=np.where(v>=0,'royalblue','darkorange'));ax.set_aspect('equal');ax.axis('off')
    s=np.interp(t,np.sort(np.arccos(mu)),v[np.argsort(np.arccos(mu))])  # sampled density vs polar angle from the bodies
    a=np.array([(2*k+1)/2*np.trapezoid(s*eval_legendre(k,np.cos(t))*np.sin(t),t) for k in range(8)]);p=a**2/(2*np.arange(8)+1);p/=p.max()
    b=fig.add_axes([0.31+c*0.5,0.785-r*0.24,0.15,0.105]);b.bar(range(8),p,color='tab:blue')
    b.set_xticks(range(8));b.set_yticks([0,1]);b.set_ylim(0,1.15);b.tick_params(labelsize=6);b.set_xlabel('Multipole L',fontsize=6)
fig.savefig(__file__.replace('.py','.png'),dpi=110,facecolor='white')
