import json,numpy as np,csv
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':18})
from matplotlib.lines import Line2D
c=json.load(open('curves.json'));g=np.array(c['grid'])
fig,axs=plt.subplots(2,1,figsize=(16,14),facecolor='white',sharex=True)
for ax,k,m,lab in ((axs[0],'A','AF','A/A_all'),(axs[1],'B','BL','B/B_all')):
  ax.axvspan(15,80,color='gold',alpha=.18);ax.axvspan(260,335,color='lightsteelblue',alpha=.35)
  ax.plot(g,c[m],color='purple',lw=3)
  ax.plot(g,c[k+'0'],color='k',lw=3)
  for kk,col in ((3,'tab:blue'),(5,'tab:green'),(10,'tab:red')):ax.plot(g,c[('YA' if k=='A' else 'PB')+str(kk)],color=col,lw=2.3,ls='--')
  ax.set_ylabel(lab);ax.set_xlim(0,360)
axs[0].set_title('A/A_all: top-k removed, ranked by joint influence on A excess + B dip over ℓ=15–80°')
axs[1].set_title('B/B_all: top-k removed, ranked by influence on mean B deficit over ℓ=260–335°')
axs[1].set_xlabel('Galactic longitude ℓ [deg]')
h=[Line2D([],[],color='k',lw=3,label='all 313 clusters'),Line2D([],[],color='tab:blue',lw=2.3,ls='--',label='top 3 removed'),Line2D([],[],color='tab:green',lw=2.3,ls='--',label='top 5 removed'),Line2D([],[],color='tab:red',lw=2.3,ls='--',label='top 10 removed'),
 Line2D([],[],color='purple',lw=3,label='model (A: AF_smooth, B: sheet BL)'),plt.Rectangle((0,0),1,1,color='gold',alpha=.3,label='ℓ≈50° deviation (ℓ=15–80°)'),plt.Rectangle((0,0),1,1,color='lightsteelblue',alpha=.5,label='ℓ≈305° deviation (ℓ=260–335°)')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False)
fig.tight_layout(rect=(0,0.09,1,1));fig.savefig('l50_jackknife_topk.png',dpi=75,facecolor='white')
