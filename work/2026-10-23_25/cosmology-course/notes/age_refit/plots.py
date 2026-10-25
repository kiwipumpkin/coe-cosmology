import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.integrate import quad
o=json.load(open('refit.json'));c=299792.458;Mpc=3.0857e22;G=6.674e-11
def phys(d,H):
    D=d['v_o']*c/H; d1=D/d['s']; M=d['mu']*(c*1e3)**2*d1*Mpc/G
    return D,M,d['ageH']*977.8/H
rows=[('baseline, H=70',o['baseline'],70),('baseline geometry, H floated',o['baseline'],o['baseline']['H_needed_for_13']),('refit, H=70 fixed, age>=13',o['constrained_H70_age13'],70)]
tab=[["Quantity"]+[r[0] for r in rows]]
P=[phys(d,H) for _,d,H in rows]
def L(n,f): tab.append([n]+[f(d,H,p) for (_,d,H),p in zip(rows,P)])
L('H at us (km/s/Mpc)',lambda d,H,p:f"{H:.1f}")
L('Observer-frame age (Gyr)',lambda d,H,p:f"{p[2]:.2f}")
L('D now (Gpc)',lambda d,H,p:f"{p[0]/1e3:.2f}")
L('Central mass M (kg)',lambda d,H,p:f"{p[1]:.2e}")
L('mu = GM/(c^2 d1)',lambda d,H,p:f"{d['mu']:.4f}")
L('beta1',lambda d,H,p:f"{d['beta1']:.3f}")
L('q',lambda d,H,p:f"{d['q']:.3f}")
L('v_o now / c',lambda d,H,p:f"{d['v_o']:.3f}")
L('CMB cancellation %',lambda d,H,p:f"{d['cancel_pct']:.1f}")
L('Residual RMS (x 1e-5 floor)',lambda d,H,p:f"{d['floor_units']:.0f}")
L('Migkas sigma / R2',lambda d,H,p:"same: 5.3-5.5σ, R² 0.743")
json.dump([[r[0],r[2]]+list(p) for r,p in zip(rows,P)],open('phys.json','w'))
import textwrap
fig,ax=plt.subplots(figsize=(16,9));ax.axis('off')
t=ax.table(cellText=[[textwrap.fill(x,24) for x in r] for r in tab],loc='center',cellLoc='center')
t.auto_set_font_size(False);t.set_fontsize(15);t.scale(1,3.2)
ax.set_title('CoE age refit: before vs after (observer-frame proper time)',fontsize=22)
plt.savefig('age_refit_table.png',dpi=110,bbox_inches='tight')
fig,ax=plt.subplots(figsize=(16,9))
for (n,d,H),p in zip(rows,P):
    s,mu,vo=d['s'],d['mu'],d['v_o'];U2=vo**2/(1-vo**2)
    r=np.linspace(1e-6,s,400);f=lambda x:np.sqrt(x)/np.sqrt(U2*x+2*mu*(1-x/s))
    tau=np.array([quad(f,0,x)[0] for x in r]);tau=tau/tau[-1]*p[2]
    ax.plot(tau,r/s,lw=3,label=f"{n}: {p[2]:.1f} Gyr")
ax.axvline(13,ls='--',c='k',label='13 Gyr (oldest stars)')
ax.set_xlabel('Observer proper time since CoE event (Gyr)',fontsize=18);ax.set_ylabel('a = r/r_now',fontsize=18)
ax.tick_params(labelsize=15);ax.legend(fontsize=15);ax.grid(alpha=.3)
plt.savefig('scale_factor_observer_age.png',dpi=110,bbox_inches='tight')
