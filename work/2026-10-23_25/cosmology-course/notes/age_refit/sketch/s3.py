import matplotlib.pyplot as plt, numpy as np
D=1.0;Rrec=D-0.35
fig,axs=plt.subplots(2,3,figsize=(18,12))
radii=[0.05,0.15,0.25,0.35,0.45,0.55]
angs=np.radians([0,45,90,135,180,225,270,315])
th=np.linspace(0,2*np.pi,400)
for ax,r,i in zip(axs.flat,radii,range(6)):
    ax.set_aspect('equal');ax.axis('off')
    ax.add_patch(plt.Circle((0,0),Rrec,color='#ffd6d6',zorder=0))
    ax.add_patch(plt.Circle((0,0),Rrec,fill=False,color='crimson',lw=2,ls='--'))
    ax.add_patch(plt.Circle((0,0),D,fill=False,color='gray',lw=1.5))
    ax.plot(0,0,'kx',ms=12,mew=2);ax.plot(D,0,'b*',ms=18)
    x,y=D+r*np.cos(th),r*np.sin(th);m=np.hypot(x,y)>=Rrec
    ax.plot(np.where(m,x,np.nan),np.where(m,y,np.nan),color='orange',lw=3)
    ax.plot(np.where(~m,x,np.nan),np.where(~m,y,np.nan),color='orange',lw=1.5,ls=':')
    for a in angs:
        ex,ey=D+r*np.cos(a),r*np.sin(a)
        if np.hypot(ex,ey)>=Rrec:
            ax.annotate('',xy=(D+0.02*np.cos(a),0.02*np.sin(a)),xytext=(ex,ey),arrowprops=dict(arrowstyle='->',color='green',lw=2))
        else:
            ax.plot(ex,ey,'o',color='crimson',ms=7)
    ax.set_title(f'Panel {i+1}: emission circle r = c(t-t*) = {r:.2f} D',fontsize=16)
    ax.set_xlim(-0.2,1.7);ax.set_ylim(-0.75,0.75)
fig.suptitle('SCHEMATIC: emission circles around us sending light in (green = from recombined matter; red dot = point still in plasma, no light)\nGray = shell around CoE (x), red dashed = recombination surface',fontsize=17)
plt.tight_layout(rect=[0,0,1,0.92]);plt.savefig('six.png',dpi=100)
