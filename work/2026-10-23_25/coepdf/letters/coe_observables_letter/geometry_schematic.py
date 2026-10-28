# Simple schematic: CoE, observer offset by D, observable universe as a circle around the observer.
import os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__))
fig,ax=plt.subplots(figsize=(10,4.2))
D,R=6.0,1.6
ax.add_patch(plt.Circle((D,0),R,fc='#cfe3ff',ec='b',lw=2))
ax.text(D,R+0.25,'our observable universe (radius R)',ha='center',color='b',fontsize=12)
ax.plot(0,0,'o',color='orange',ms=16,mec='k');ax.text(0,-0.6,'centre of expansion (CoE)',ha='center',fontsize=12)
for r in (1,2,3.2):ax.add_patch(plt.Circle((0,0),r,fill=False,ec='#b35900',ls='--',lw=1.8))
ax.plot(D,0,'k*',ms=18);ax.text(D,-0.45,'us',ha='center',fontsize=12)
ax.annotate('',(0.25,0),(D-0.25,0),arrowprops=dict(arrowstyle='->',lw=1.8))
ax.text(D/2,0.18,'line to the CoE, distance D  (Galactic $l\\approx180^\\circ$)',ha='center',fontsize=12)
ax.plot([D,D+R],[0,0],'b--',lw=1);ax.text(D+R/2,-0.3,'R',color='b',ha='center',fontsize=12)
ax.set_aspect('equal');ax.set_xlim(-3.5,D+R+0.5);ax.set_ylim(-2.1,2.3);ax.axis('off')
fig.tight_layout();fig.savefig(H+'/geometry_schematic.png',dpi=150,facecolor='white')
