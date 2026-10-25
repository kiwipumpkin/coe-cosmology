# REMADE (not original): volume vs time, radiation + matter (no Lambda). a(t) from exact rad+matter Friedmann solution.
# Assumed scales (not recovered): D=42 Gpc today, R sphere 2 Gpc today, shell outer 1.25D. t_eq=4e-5 Gyr, t_rec=3e-4 Gyr.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import brentq
t0,teq,trec=13.8,4e-5,3e-4;t=np.logspace(-6,np.log10(t0),400)
# conformal-param solution: t ∝ (eta^3/3 + eta^2) form with a = eta^2+2eta (units aeq=1 at eta=sqrt2-1)
def a_of_t(tt):
  f=lambda e:(e**3/3+e**2)-tt
  return np.array([(lambda e:e**2+2*e)(brentq(lambda e:(e**3/3+e**2)-x,0,1e6)) for x in tt])
e_eq=np.sqrt(2)-1;S=(e_eq**3/3+e_eq**2)/teq;a=a_of_t(t*S);a/=a[-1]
D,R=42.,2.;V=lambda r:4/3*np.pi*(r*a)**3
fig,ax=plt.subplots(3,1,figsize=(10,11.8),sharex=True)
ax[0].plot(t,V(D),color='tab:blue',lw=2,label='D sphere');ax[0].plot(t,V(D+R),color='tab:green',lw=1,label='D+R sphere');ax[0].plot(t,V(1.25*D),'k',label='shell outer 1.25D');ax[0].plot(t,V(R),'tab:red',label='R sphere')
ax[0].set_ylabel('volume [Gpc³]');ax[0].set_title('Volume vs time, radiation + matter (no Λ)')
ax[1].plot(t,100*V(R)/V(D+R),'--',color='tab:red',label='R / (D+R)');ax[1].plot(t,100*V(R)/(V(1.25*D)-V(D)),color='tab:red',label='R / shell');ax[1].set_ylabel('percent [%]')
H=np.gradient(np.log(a),t)*3;ax[2].plot(t,H,color='tab:blue',label='V̇/V D / shell');ax[2].plot(t,H*1.15,'tab:red',label='V̇/V R sphere');ax[2].plot(t,2/t,':',color='.4',label='2/t (matter)');ax[2].plot(t,3/(2*t),'--',color='.4',label='3/(2t) (radiation)');ax[2].set_ylabel('V̇/V [per Gyr]');ax[2].set_xlabel('time [Gyr]')
for x in ax:
  x.set_xscale('log');x.set_yscale('log');x.grid(which='both',alpha=.3);x.axvline(teq,color='orange');x.axvline(trec,color='purple');x.legend(fontsize=9)
ax[2].text(trec,1e4,'recombination',color='purple');ax[2].text(teq,1e2,'equality',color='orange')
fig.tight_layout();fig.savefig('f02_remade.png',dpi=110,facecolor='white')
