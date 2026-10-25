import matplotlib.pyplot as plt, numpy as np
Rrec=0.65;U=1.5;dr=0.07
fig,axs=plt.subplots(2,3,figsize=(18,12))
angs=np.radians([0,45,90,135,180,225,270,315])
for i,ax in enumerate(axs.flat):
    ax.set_aspect('equal');ax.axis('off')
    ax.add_patch(plt.Circle((0,0),Rrec,color='#ffd6d6',zorder=0))
    ax.add_patch(plt.Circle((0,0),Rrec,fill=False,color='crimson',lw=2,ls='--'))
    ax.plot(0,0,'kx',ms=12,mew=2);ax.text(0.02,-0.08,'CoE',fontsize=13)
    xu=Rrec+U*dr*i
    ax.plot(xu,0,'b*',ms=20,zorder=5)
    if i>0: ax.add_patch(plt.Circle((xu,0),dr*i,fill=False,color='blue',lw=2.5))
    if i>0:
        xp=Rrec+U*dr*(i-1); r=dr*i
        ax.plot(xp,0,'*',color='lightblue',ms=16)
        ax.text(xp-.04,.04,f'us at step {i-1}',fontsize=11,color='steelblue')
        ax.add_patch(plt.Circle((xp,0),r,fill=False,color='orange',lw=3,ls='--'))
        for a in angs:
            ex,ey=xp+r*np.cos(a),r*np.sin(a)
            ax.plot(ex,ey,'o',color='orange',ms=7)
            ax.annotate('',xy=(xu+0.02*np.cos(np.arctan2(ey,ex-xu)),0.02*np.sin(np.arctan2(ey,ex-xu))),xytext=(ex,ey),arrowprops=dict(arrowstyle='->',color='green',lw=1.8))
    ax.set_title(f'Step {i}: us at {xu:.2f} from CoE\nseeing light emitted at step {max(i-1,0)}',fontsize=15)
    ax.set_xlim(-0.2,1.9);ax.set_ylim(-0.6,0.6)
fig.suptitle('SCHEMATIC: light seen at each step (green) was emitted one step earlier from the circle\naround where we were then (orange dashed, light-blue star); blue solid = current-step circle around us now (U=1.5c)',fontsize=17)
plt.tight_layout(rect=[0,0,1,0.9]);plt.savefig('six_v5.png',dpi=100)
