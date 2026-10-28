# Galactic Mollweide map of the published dipoles/axes and anomaly directions, with the CoE direction.
import os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__))
P=[('CMB dipole',264.0,48.3,'r','o',True),('CatWISE dipole',238.2,28.8,'m','s',True),('NVSS dipole',196,46,'purple','D',True),
('Bulk flow',287,8,'g','^',True),('Axis of Evil',260,60,'k','*',True),
('Migkas max. anisotropy',280,-15,'orange','P',False),('CMB Cold Spot',209,-57,'b','h',False),('Loop I / NPS',329,17.5,'brown','X',False),('CoE (model)',180,-15,'gold','o',True)]
def xy(l,b):return -np.radians(((l+180)%360)-180),np.radians(b)
fig=plt.figure(figsize=(12,6.6));ax=fig.add_subplot(111,projection='mollweide');ax.grid(alpha=.4)
for n,l,b,c,mk,anti in P:
  x,y=xy(l,b);ax.plot(x,y,mk,color=c,mec='k',ms=14 if mk!='*' else 20,label=f'{n} ({l:g}°, {b:+g}°)')
  if anti:x2,y2=xy(l+180,-b);ax.plot([x,x2],[y,y2],'-',color=c,lw=2,alpha=.8);ax.plot(x2,y2,mk,mfc='none',mec=c,mew=2,ms=12 if mk!='*' else 18)
ax.set_xticks(np.radians([-150,-120,-90,-60,-30,0,30,60,90,120,150]));ax.set_xticklabels(['150','120','90','60','30','0','330','300','270','240','210'])
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.05),ncol=3,fontsize=9,frameon=False)
ax.set_title('Galactic Mollweide, l=0 at centre, increasing to the left. Filled: published direction; open: antipode; line joins each pair',fontsize=10,pad=14)
fig.tight_layout();fig.savefig(H+'/dipoles_sky_map_lines.png',dpi=140,facecolor='white',bbox_inches='tight')
