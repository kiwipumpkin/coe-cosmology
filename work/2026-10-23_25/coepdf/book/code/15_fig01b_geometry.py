# Fig 1 bottom: observer O (Earth) inside thick shell, distance D from CoE core, vision sphere R.
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,numpy as np
from matplotlib.patches import Circle,FancyArrowPatch
plt.rcParams.update({'mathtext.fontset':'cm','font.family':'serif','font.size':15})
fig=plt.figure(figsize=(5.34,4.12),dpi=200,facecolor='white');ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,481);ax.set_ylim(371,0);ax.axis('off')
C=(38,185);O=(214,185);D=O[0]-C[0]
ax.add_patch(Circle(C,262,fc='lightgrey',ec='0.35',lw=1.5))
ax.add_patch(Circle(C,96,fc='white',ec='0.35',lw=1.5))
ax.add_patch(Circle(C,D,fill=False,ec='tab:blue',ls='--',lw=2))
ax.add_patch(Circle(O,76,fc='tab:green',alpha=.3,ec='tab:green',lw=1.5))
ar=lambda a,b,c,lw=1.6:ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=14,color=c,lw=lw))
# D: CoE -> observer, drawn as chord-free radius at an angle ending on dashed circle (|D| = CoE-to-O distance)
ar((C[0]+6,C[1]),(O[0]-7,O[1]),'tab:blue',2);ax.text(118,C[1]-8,'$D$',color='tab:blue',fontsize=18,ha='center')
t2=np.radians(125);ar(O,(O[0]+76*np.cos(t2),O[1]-76*np.sin(t2)),'tab:green',2)
ax.text(O[0]-42,O[1]-20,'$R$',color='tab:green',fontsize=17)
ar(O,(472,185),'#0055cc',2.5);ax.text(390,218,r'$V_O\!\approx\!7.7c$',color='#0055cc',fontsize=15)
ax.plot(*C,'o',color='tab:red',ms=10);ax.plot(*O,'o',color='#0055cc',mec='k',ms=9)
ax.text(8,C[1]+30,'CoE core',color='tab:red',ha='left',fontsize=15)
ax.text(68,118,'empty\ninterior',color='0.4',ha='center',fontsize=14,linespacing=1.1)
ax.text(O[0]+6,O[1]+40,'Earth\n(observer, $O$)',color='k',fontsize=13,ha='left',linespacing=1.05,bbox=dict(fc='white',ec='none',alpha=.75,pad=1))
ax.text(200,40,'shell',rotation=-42,color='0.3',fontsize=15,ha='center',va='center')
ax.text(270,40,'core $r_s$ = 163 Gpc\n(4.9 D), off scale',color='tab:red',fontsize=14)
ax.add_patch(plt.Rectangle((0,0),481,371,fill=False,ec='k',lw=2))
fig.savefig('fig/Fig01b_geometry.png',facecolor='white')
