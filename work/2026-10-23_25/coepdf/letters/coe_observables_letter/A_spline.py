# Periodic cubic spline through Migkas A/A_all central values (and max/min envelopes) vs Galactic longitude.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high']
o=np.argsort(l);l,y,el,eh=l[o],y[o],el[o],eh[o]
def sp(v):return CubicSpline(np.r_[l,l[0]+360],np.r_[v,v[0]],bc_type='periodic')
g=np.linspace(0,360,36001)
fig,ax=plt.subplots(figsize=(11,5));ax.errorbar(l,y,[el,eh],fmt='ko',ms=3,capsize=2,label='Migkas et al. (2020) $A/A_{\\rm all}$')
ax.plot(g,sp(y)(g%360+0*g) if False else sp(y)(np.where(g<l[0],g+360,g)),'b-',label='periodic cubic spline, central')
for v,c,n in [(y+eh,'crimson','max envelope'),(y-el,'royalblue','min envelope')]:ax.plot(g,sp(v)(np.where(g<l[0],g+360,g)),'--',color=c,lw=1,label='spline, '+n)
ax.set(xlabel='Galactic longitude l [deg]',ylabel='$A/A_{\\rm all}$',xlim=(0,360));ax.set_xticks(range(0,361,60));ax.grid(alpha=.3);ax.legend(fontsize=8)
plt.tight_layout();plt.savefig('A_spline.png',dpi=140)
print(l[0],l[-1],len(l))
