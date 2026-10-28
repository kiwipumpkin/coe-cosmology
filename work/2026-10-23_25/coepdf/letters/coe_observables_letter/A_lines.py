# Migkas A/A_all: straight lines connecting the data points (central, max, min), periodic wrap at 360.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high']
o=np.argsort(l);l,y,el,eh=l[o],y[o],el[o],eh[o]
g=np.linspace(0,360,36001);f=lambda v:np.interp(g,l,v,period=360)
fig,ax=plt.subplots(figsize=(11,5));ax.errorbar(l,y,[el,eh],fmt='ko',ms=3,capsize=2,label='Migkas et al. (2020) $A/A_{\\rm all}$')
ax.plot(g,f(y),'b-',label='connect the dots, central')
ax.plot(g,f(y+eh),'--',color='crimson',lw=1,label='connect the dots, max');ax.plot(g,f(y-el),'--',color='royalblue',lw=1,label='connect the dots, min')
ax.set(xlabel='Galactic longitude l [deg]',ylabel='$A/A_{\\rm all}$',xlim=(0,360));ax.set_xticks(range(0,361,60));ax.grid(alpha=.3);ax.legend(fontsize=8)
plt.tight_layout();plt.savefig('A_lines.png',dpi=140)
