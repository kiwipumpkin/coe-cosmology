import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
exec(open('/tmp/spots.py').read());M=json.load(open('/tmp/M25.json'));l=np.array(M['l'])
Z={'A':[22,99,262,310],'B':[14,104,266,321]}
fig,axs=plt.subplots(1,2,figsize=(20,11),subplot_kw={'projection':'polar'})
for ax,k in zip(axs,'AB'):
  ax.set_theta_zero_location('N');ax.set_theta_direction(1)
  L=np.linspace(0,360,721);v=np.interp(L,l,M[k],period=360);V=np.abs(v).max()
  th=np.radians(L)
  for i in range(len(L)-1):
    c=plt.cm.bwr(.5+.5*v[i]/V)
    ax.fill_between(th[i:i+2],0.55,1.0,color=c,lw=0)
  for z in Z[k]: ax.plot([np.radians(z)]*2,[0.5,1.05],'k-',lw=2)
  for n,(nm,sl,sb,cc) in enumerate(SPOTS[k]):
    val=float(np.interp(sl,l,M[k],period=360))
    zone=None
    zs=sorted(Z[k]);
    ok=(val>0)==(cc=='hot')
    t=np.radians(sl);col=CC[cc]
    ax.plot([t],[0.78],'o',ms=24,mfc=col,mec='k',mew=2,zorder=5)
    ax.text(t,0.78,str(n+1),ha='center',va='center',color='w',weight='bold',fontsize=13,zorder=6)
    ax.text(t,1.22,f'{nm}\nℓ={sl}°  {val:+.1f}%',ha='center',va='center',fontsize=12,color=col,weight='bold',
      bbox=dict(fc='w',ec=col,lw=1.5))
  ax.set_ylim(0,1.35);ax.set_yticks([]);ax.set_xticks(np.radians(range(0,360,30)))
  ax.set_xticklabels([f'{d}°' for d in range(0,360,30)],fontsize=12);ax.grid(alpha=.3)
  ax.set_title(f'{k} residual vs ℓ (±25° mean)  —  ring colour: red > 0, blue < 0; black = zero crossings',fontsize=15,pad=45)
  sm=plt.cm.ScalarMappable(cmap='bwr',norm=plt.Normalize(-V,V))
  cb=fig.colorbar(sm,ax=ax,orientation='horizontal',fraction=.04,pad=.08);cb.set_label(f'{k} residual (%)',fontsize=12)
fig.legend([Line2D([],[],marker='o',ls='',mfc='red',mec='k',ms=14),Line2D([],[],marker='o',ls='',mfc='blue',mec='k',ms=14)],
  ['hot spot (needs > 0)','cold spot (needs < 0)'],loc='lower center',ncol=2,fontsize=13,frameon=False)
fig.tight_layout(rect=(0,.05,1,1));fig.savefig('polar_AB.png',dpi=110,facecolor='white')
